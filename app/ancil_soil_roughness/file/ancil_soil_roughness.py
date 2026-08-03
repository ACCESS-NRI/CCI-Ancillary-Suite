#!/usr/bin/env python
"""
Soil roughness application
**************************

Derive the soil roughness ancillary by:

- Regrid the provided LAI source climatology to the source roughness grid.
    - The soil roughness is dependent on the LAI at source resolution.
- Masking the roughness source by the following conditions:
    - Roughness greater than 0.08 AND where the annual mean LAI is > 4
    - The annual minimum LAI is > 2
- Regrid the soil roughness source to the provided LAI grid by assuming a
  blending height of 10m (bespoke regrid).
- Data points which are not available in the raw source data (as opposed to
  the LAI filter) are filled in using a spiral search interpolation.

Options:
 * --target-lsm : The land/sea mask of the target grid
 * --leaf-area-index : A leaf area index file. This should be the intermediate
   pre-processed file for an LAI climatology, so it is NOT on the target 
   grid, and is NOT on land surface tiles
 * --outfile : the output file
Arguments:
 * The source data file follows the options as a single argument.

Example Usage:
ancil_soil_roughness.py --target-lsm ${lsm} --leaf-area-index ${lai} \
                        --output ${outfile} ${srcfile}

Application contact: Melissa Brooks (Melissa.Brooks@metoffice.gov.uk)

"""

# TODO:
#
# - Revisit hardcoding this blending height (10m) as part of the source roughness
#   source dataset and interpret it here.
# - Can the bespoke regrid strategy be replaced with a centralised one, or a
#   change of workflow allow the use of existing methods?


import ants
import ants.decomposition as decomp
import ants.io.save as save
import iris
import numpy as np

iris.FUTURE.save_split_attrs = True

def main(source_path, lai_path, target_lsm, outfile, netcdf_only):
    "main top/level"
    z0_src, lsm, lai = load(source_path, target_lsm, lai_path)
    res = decomp.decompose(gen_soil_roughness, [z0_src, lai], lsm)

    save.netcdf(res, outfile)
    if not netcdf_only:
        save.ancil(res, outfile)


def load(source_path, lsm_path, lai_path):
    """
    loads in the required source data, the target land/sea mask and
    the leaf area index
    """
    lai_path = ants.config.filepath_readable(lai_path)

    z0_constr = iris.AttributeConstraint(STASH="m01s00i097")
    z0_src = ants.io.load.load_cube(source_path, z0_constr)

    lsm = ants.io.load.load_cube(lsm_path, "land_binary_mask")
    lsm = lsm.copy(lsm.data.astype("bool", copy=False))

    lai = ants.io.load.load_cube(lai_path)

    return z0_src, lsm, lai


def gen_soil_roughness(source, lsm):
    """
    produces soil roughness from the source roughness lenght and lai data
    to the required lsm
    """
    z0_src, lai = source
    # now we want to produce the time-mean and time-minimum LAI, which are then
    # regridded to the source resolution
    lai_tmean, lai_tmin = time_mean_min_and_regrid(lai, z0_src)

    # fill value for z0, if just filling in a background value
    # (1e-3 was the GA6 value):
    z0_fill_value = 1.00000e-3
    blending_height = 10

    # going to regrid the data in three stages (once raw,
    # once using the masked data to get 'known good' data and finally 
    # using that with the masked data filled in.
    z0_raw = regrid_z0_cubes(z0_src, lsm, blending_height)

    mask = gen_z0_mask(z0_src, lai_tmean, lai_tmin)
    masked_z0_src = z0_src.copy()
    masked_z0_src.data[mask] = np.ma.masked
    z0_below_lai = regrid_z0_cubes(masked_z0_src, lsm, blending_height)

    # this is the more stringent application for good data...
    known_good_data_ind = np.where(np.logical_not(z0_below_lai.data.mask))
    known_good_data = z0_below_lai.data[known_good_data_ind]

    # now do the same as before but this time fill in any missing data first.
    masked_z0_src.data[mask] = z0_fill_value
    z0_filled_bl_lai = regrid_z0_cubes(masked_z0_src, lsm, blending_height)

    # Now use the filled cube as a base, and overwrite with the points that
    # had enough obs:

    # if lsm says there is ocean, set as z0 as masked
    z0_filled_bl_lai.data[~lsm.data] = np.ma.masked
    # if lsm says there is land but z0 is masked, set to our fill value
    z0_filled_bl_lai.data[
        lsm.data & np.ma.getmaskarray(z0_filled_bl_lai.data)
    ] = z0_fill_value

    # also overwrite with the 'known_good_data' where possible:
    z0_filled_bl_lai.data[known_good_data_ind] = known_good_data

    # now find points that had no data at all, so it's all fill value:
    all_fill = z0_filled_bl_lai.data.data == z0_fill_value
    # if they're all fill because the raw data was missing for those points
    # then we need to do an interpolation to fill those in:
    missing_from_raw = np.where(np.logical_and(all_fill, z0_raw.data.mask))
    # now those points are missing, to be filled:
    if np.size(missing_from_raw) > 0:
        # mask those points out:
        z0_filled_bl_lai.data.mask[missing_from_raw] = True
        is_land = True # flag for make_consistent_with_lsm
        ants.analysis.make_consistent_with_lsm([z0_filled_bl_lai], lsm, is_land)

    out = z0_filled_bl_lai

    return out


def time_mean_min_and_regrid(lai, z0_src):
    """
    Prepares the input LAI by producing a time mean and min, and regrids
    to the z0_src grid.

    Unsegmented, this peaks at 8GB RAM
    """

    lai_tmean_sgrid = lai.collapsed("time", iris.analysis.MEAN, mdtol=1.0)
    lai_tmean = lai_tmean_sgrid.regrid(z0_src, iris.analysis.Linear())

    lai_tmin_sgrid = lai.collapsed("time", iris.analysis.MIN, mdtol=1.0)
    lai_tmin = lai_tmin_sgrid.regrid(z0_src, iris.analysis.Linear())
    return lai_tmean, lai_tmin


def gen_z0_mask(z0_src, lai_tmean, lai_tmin):
    """
    Generate mask for z0_src and regrids a z0 cube. Returns True if the source
    z0 data is not to be used.

    """
    # This is what high z0 is defined as (m):
    high_z0_value = 0.08
    mask = z0_src.data > high_z0_value

    # Threshold any lai level (time mean) above which Prigent data is ignored.
    lai_mean_t = 4.0
    mask &= lai_tmean.data > lai_mean_t

    # Threshold any lai level (time min) above which Prigent data is ignored.
    lai_min_t = 2.0
    mask |= lai_tmin.data > lai_min_t
    return mask


def regrid_z0_cubes(z0_in_cube, grid_cube, blending_height):
    """Appropriate regridding of z0, from one IRIS cube to another.
    This uses area weighted averaging of the flux weighted z0
    (assuming a blending height)
    """
    # convert to fz0
    fz0_cube = _z0_cube_to_fz0(z0_in_cube, blending_height)

    # then regrid, and convert back to z0
    scheme = ants.regrid.rectilinear.TwoStage(mdtol=0.9)
    scheme = ants.regrid.GeneralRegridScheme(horizontal_scheme=scheme)
    rgd = fz0_cube.regrid(grid_cube, scheme)
    z0_out = _fz0_cube_to_z0(rgd, blending_height)

    # convert the units back, if it's been changed in the conversion:
    if z0_out.units != z0_in_cube.units:
        z0_out.convert_units(z0_in_cube.units)

    return z0_out


def _z0_cube_to_fz0(z0_in_cube, blending_height):
    "converts a z0 cube to a flux weighted z0 (assuming a blending height)"

    if z0_in_cube.units != "m":
        use_z0 = z0_in_cube.copy()
        use_z0.convert_units("m")
    else:
        use_z0 = z0_in_cube
    # first convert input data to fz0
    fz0_cube = use_z0.copy()
    fz0_cube.data = 1.0 / (np.log(blending_height / fz0_cube.data)) ** 2.0

    return fz0_cube


def _fz0_cube_to_z0(fz0_in_cube, blending_height):
    "converts a flux weighted z0 cube to a z0 (assuming a blending height)"
    z0_cube = fz0_in_cube.copy()
    # then convert back:
    z0_cube.data = blending_height * np.exp(-np.sqrt(1.0 / z0_cube.data))
    return z0_cube


if __name__ == "__main__":
    parser = ants.AntsArgParser(target_lsm=True)
    parser.add_argument(
        "--leaf-area-index",
        type=str,
        help="Path to leaf area index file",
        required=True,
    )
    args = parser.parse_args()
    main(args.sources, args.leaf_area_index, args.target_lsm, args.output, args.netcdf_only)
