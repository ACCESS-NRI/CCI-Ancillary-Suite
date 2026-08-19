# The first functional version of this script was developed by hand. The rework
# to take advantage of dask parallelism was written by Claude Sonnet 5, with
# the original script as input and direction to make it parallelise using dask.

import functools
import math

import xarray
import numpy
import argparse
import scipy.ndimage
import dask.array as da
from dask.distributed import Client
from dask.diagnostics import ProgressBar


def _parse_args():
    parser = argparse.ArgumentParser(
            """Script to merge the CCI land cover for a specific year with the
            permanent water bodies dataset. The land cover dataset is at 300m
            resolution, while the permanent water is at 150m resolution."""
            )
    parser.add_argument(
            '--CCI-land-cover',
            type=str,
            help='Source CCI land cover dataset',
            required=True
            )
    parser.add_argument(
            '--CCI-water-dir',
            type=str,
            help='Directory containing the CCI water bodies sources',
            required=True
            )
    parser.add_argument(
            '-o',
            '--output',
            help='Output path',
            required=True
            )
    parser.add_argument(
            '--nchunks',
            type=int,
            default=16,
            help=(
                'Target number of chunks to split the land cover grid into '
                'for parallel processing (roughly one chunk per available '
                'worker thread is a reasonable starting point). The grid is '
                'split along latitude and longitude in proportion to its '
                'aspect ratio to keep chunks roughly square; the actual '
                'chunk count may differ slightly from this value due to '
                'rounding. The water bodies datasets are chunked at twice '
                'the resulting cell size along each dimension so that chunk '
                'boundaries line up exactly with the 2x2 upscaling.'
                )
            )
    parser.add_argument(
            '--fill-depth',
            type=int,
            default=64,
            help=(
                'Maximum search radius, in permanent-water (150m) cells, used '
                'when filling land-cover cells that disagree with the '
                'permanent water mask. This bounds the nearest-neighbour fill '
                'to a halo around each chunk so it can run in parallel; it '
                'must be larger than the largest expected gap in your data or '
                'cells near the middle of a large misclassified patch will be '
                'filled from the wrong side of the chunk.'
                )
            )
    parser.add_argument(
            '--n-workers',
            type=int,
            default=4,
            help='Number of Dask worker processes'
            )
    parser.add_argument(
            '--threads-per-worker',
            type=int,
            default=2,
            help='Threads per Dask worker'
            )

    return parser.parse_args()


def main(
        land_cover_source,
        water_bodies_dir,
        output,
        nchunks,
        fill_depth,
        n_workers,
        threads_per_worker,
        ):
    client = Client(n_workers=n_workers, threads_per_worker=threads_per_worker)

    # Load in the data sources. These aren't chunked yet - the chunking is
    # done explicitly below so we can guarantee the water bodies chunks are
    # exactly double the land cover chunks, and aligned to the same origin.
    land_cover = xarray.open_dataset(land_cover_source)

    # The water bodies have 2 datasets which have fixed names
    ocean_water = xarray.open_dataset(water_bodies_dir + 'ESACCI-LC-L4-WB-Ocean-Map-150m-P13Y-2000-v4.0.tif')
    all_water = xarray.open_dataset(water_bodies_dir + 'ESACCI-LC-L4-WB-Map-150m-P13Y-2000-v4.0.tif')

    land_cover, ocean_water, all_water = _rechunk_matching(
            land_cover, ocean_water, all_water, nchunks
            )

    upscaled_water_bodies = upscale_water_bodies(ocean_water, all_water)

    land_cover = overlay_water_on_land_cover(
            upscaled_water_bodies, land_cover, fill_depth
            )

    land_cover = drop_unnecessary_data(land_cover)
    add_coord_system(land_cover)
    correct_metadata(land_cover)

    # to_netcdf will trigger computation of the whole (still lazy) dask graph.
    # compute=False + an explicit .compute() lets us wrap it in a progress bar.
    delayed_write = land_cover.to_netcdf(output, compute=False)
    with ProgressBar():
        delayed_write.compute()

    client.close()


def _compute_chunk_sizes(lat_size, lon_size, nchunks):
    """
    Given the full extent of the land cover grid and a target number of
    chunks, work out a chunk size (in cells) for each spatial dimension such
    that the grid splits into approximately `nchunks` chunks in total. The
    split between dimensions is weighted by the grid's aspect ratio so
    chunks come out roughly square, rather than being long thin strips.
    """
    nchunks = max(1, nchunks)
    aspect = lon_size / lat_size

    n_lat = max(1, round(math.sqrt(nchunks / aspect)))
    n_lon = max(1, round(nchunks / n_lat))

    lat_chunk = math.ceil(lat_size / n_lat)
    lon_chunk = math.ceil(lon_size / n_lon)

    return lat_chunk, lon_chunk


def _rechunk_matching(land_cover, ocean_water, all_water, nchunks):
    """
    Chunk the land cover dataset into approximately `nchunks` chunks, and the
    two (double resolution) water bodies datasets at exactly twice the
    resulting cell size along each dimension, so that every land cover chunk
    corresponds exactly to a 2x2 block of water bodies chunks. Dimension
    names are read from each dataset rather than assumed, since the land
    cover (netCDF) and water bodies (GeoTIFF) sources typically use different
    dimension naming conventions.
    """
    lat_dim, lon_dim = land_cover['lccs_class'].dims[-2:]
    lat_size = land_cover.sizes[lat_dim]
    lon_size = land_cover.sizes[lon_dim]

    lat_chunk, lon_chunk = _compute_chunk_sizes(lat_size, lon_size, nchunks)

    land_cover = land_cover.chunk({lat_dim: lat_chunk, lon_dim: lon_chunk})

    y_dim, x_dim = ocean_water['band_data'].dims[-2:]
    ocean_water = ocean_water.chunk({y_dim: lat_chunk * 2, x_dim: lon_chunk * 2})
    all_water = all_water.chunk({y_dim: lat_chunk * 2, x_dim: lon_chunk * 2})

    return land_cover, ocean_water, all_water


def upscale_water_bodies(ocean_water, all_water):
    """
    Combine and upscale the provided water bodies dataset to create a single
    dataset at the desired 300m resolution, which differentiates between ocean
    and inland water. Returned dataset has 0 for ocean, 1 for land, 2 for
    inland water.

    This is embarrassingly parallel over chunks: each 300m output chunk
    depends only on the corresponding 2x2 block of 150m input chunks, which
    is guaranteed by `_rechunk_matching`. All mutation is expressed with
    `dask.array.where` instead of boolean-mask assignment, since dask arrays
    don't support efficient in-place item assignment the way numpy arrays do.
    """

    # First, combine the datasets to differentiate between ocean and inland
    # water. The ocean dataset has 0 = ocean, 1 = other, while the all water
    # dataset has 1 = other, 2 = water.
    ocean_data = ocean_water["band_data"][0, :, :].data
    all_data = all_water["band_data"][0, :, :].data

    water_data = da.where(all_data == 1, 1, 0)
    water_data = da.where(
            da.logical_and(all_data == 2, ocean_data == 1), 2, water_data
            )

    # We should be able to do this somewhat cleverly, instead of relying on
    # some windowed most common algorithm. Given we know the only possible
    # values are 0, 1, 2, we apply a three step process:
    # 1. Change the value "2" to "5" everywhere
    # 2. Perform the windowed sum, so we get the sum of values in every 2x2 box
    # 3. Inspect the resulting value, to determine what the most common value
    #   was in the original 2x2 box.
    # The bias will point "downwards"; the lowest flag value will take
    # precedence when two flags are equally most common.
    # This means the rules when inspecting the resulting values are:
    # 1. If value > 10, then the most common flag was a 5 (originally 2)
    # 2. If mod(value, 5) > 2, then the most common flag was a 1
    # 3. Otherwise, the most common flag was a 0.
    water_data = da.where(water_data == 2, 5, water_data)

    upscaled_data = water_data[0::2, 0::2] + water_data[1::2, 0::2] + \
            water_data[0::2, 1::2] + water_data[1::2, 1::2]

    flag_data = da.zeros_like(upscaled_data)
    flag_data = da.where(upscaled_data > 10, 2, flag_data)
    flag_data = da.where(da.mod(upscaled_data, 5) > 2, 1, flag_data)

    return flag_data


def overlay_water_on_land_cover(water, land_cover, fill_depth):
    """
    Apply the permanent water map over the top of the land cover map. The
    permanent water map is taken as the absolute truth in regards to water, and
    will override any conflicting land values in the land cover map. Any cells
    considered water by the land cover map, but land by the permanent water
    map, will be filled by the nearest valid land point i.e. nearest point
    along the same latitude considered land by the permanent water map.

    The nearest-valid-point fill is done per-chunk with a `fill_depth`-cell
    halo via `map_overlap`, since a true global distance transform can't be
    parallelized without materializing the whole array. See `--fill-depth`
    for the trade-off this introduces.
    """

    flag_meanings = land_cover['lccs_class'].attrs['flag_meanings'].split()
    flag_values = land_cover['lccs_class'].attrs['flag_values']

    # Identify the flag associated with water
    water_index = flag_meanings.index('water')
    water_flag = flag_values[water_index]

    # Split the water into sea_ocean_water and inland_water. To make it clear
    # that this is a different flag to the original water, don't reuse the same
    # flag- use the next two free flag values
    sea_ocean_water_flag = water_flag + 1
    inland_water_flag = water_flag + 2

    # Now insert them in the locations adjacent to the original water flag
    flag_meanings.insert(water_index + 1, 'sea_ocean_water')
    flag_meanings.insert(water_index + 2, 'inland_water')
    flag_values = numpy.insert(flag_values, water_index + 1, sea_ocean_water_flag)
    flag_values = numpy.insert(flag_values, water_index + 2, inland_water_flag)

    # Remove the now-unused water flag
    flag_meanings.pop(water_index)
    flag_values = numpy.delete(flag_values, water_index)

    land_cover['lccs_class'].attrs['flag_meanings'] = ' '.join(flag_meanings)
    land_cover['lccs_class'].attrs['flag_values'] = flag_values

    # Apply the new flag values for sea_ocean_water and inland_water to the
    # land cover. Any remaining water flag values are cells that we want to
    # write over with the nearest non-water value.
    lccs = land_cover['lccs_class'].data
    lccs = da.where(water == 0, sea_ocean_water_flag, lccs)
    lccs = da.where(water == 2, inland_water_flag, lccs)

    fill_block = functools.partial(_fill_nearest_block, water_flag=water_flag)
    lccs = lccs.map_overlap(
            fill_block,
            depth=fill_depth,
            boundary='reflect',
            dtype=lccs.dtype,
            )

    land_cover['lccs_class'].data = lccs
    return land_cover


def _fill_nearest_block(block, water_flag):
    """
    Fill cells equal to `water_flag` with the value of the nearest cell that
    isn't `water_flag`, operating on a single chunk (plus its overlap halo).
    Runs scipy's distance transform locally, which is why it needs the halo
    from `map_overlap` to see past its own chunk boundary.
    """
    mask = block == water_flag
    if not mask.any():
        return block

    _, indices = scipy.ndimage.distance_transform_edt(mask, return_indices=True)
    return block[tuple(indices)]


def drop_unnecessary_data(land_cover):
    """
    The land cover dataset includes ancillary variables which are not
    necessary for the land cover mapping, and consume significant disk space.
    Drop them from the dataset.
    """
    land_cover["lccs_class"].attrs.pop("ancillary_variables")
    land_cover = land_cover[["lccs_class", "lat", "lon"]]

    return land_cover


def add_coord_system(dataset):
    """
    ANTS/Iris requires a coordinate system to correctly interpret the data.
    """
    dataset["latitude_longitude"] = xarray.DataArray(
            0,
            attrs={
                "grid_mapping_name": "latitude_longitude",
                "longitude_of_prime_meridian": 0.0,
                "earth_radius": 6371229.0
                }
            )

    dataset["lccs_class"].attrs["grid_mapping"] = "latitude_longitude"


def correct_metadata(dataset):
    """
    The dataset includes a no_data flag, but there are no instances of no_data
    in the dataset. Remove it from the flag specifications. The data is assigned
    a valid range based of unsigned integers, but the flags are based on signed
    integers. The data is  designated as _Unsigned, but the flags are provided
    as signed integers. Drop this attribute to ensure it's handled correctly.
    """
    dataset["lccs_class"].attrs["flag_values"] = \
            dataset["lccs_class"].attrs["flag_values"][1:]
    dataset["lccs_class"].attrs["flag_meanings"] = \
            " ".join(dataset["lccs_class"].attrs["flag_meanings"].split(" ")[1:])

    dataset["lccs_class"].attrs["valid_min"] = -127
    dataset["lccs_class"].attrs["valid_max"] = 128
    dataset["lccs_class"].encoding.pop("_Unsigned")


if __name__ == "__main__":
    args = _parse_args()

    main(
            args.CCI_land_cover,
            args.CCI_water_dir,
            args.output,
            args.nchunks,
            args.fill_depth,
            args.n_workers,
            args.threads_per_worker,
            )
