# (C) Crown Copyright, Met Office. All rights reserved.
"""
Common functions used for the orographic wavedrag and radiation parameter processing.

"""
import ants
import cartopy.crs as ccrs
import iris
import iris.analysis.calculus
import numpy as np
from scipy import ndimage, signal


def calc_mean(source, target):
    """
    Calculate the mean field of the orography.

    Calculates the mean field from the source orography by regridding it onto the
    target. Also fixes the metadata (name and STASH) of the field.

    Parameters
    ----------
    source : :class:`iris.cube.Cube`
        Orography source dataset.
    target : :class:`iris.cube.Cube`
        Target landsea mask.

    Returns
    -------
    : :class:`iris.cube.Cube`
        orography mean height (m01s00i033).

    """
    mean = source.regrid(
        target, ants.regrid.GeneralRegridScheme(horizontal_scheme="TwoStage")
    )
    mean.rename("orography mean height")
    mean.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i033")
    return mean


def gen_ocean_mask(landseamask):
    """
    Calculate the open ocean mask.

    Extracts the open ocean mask by flood filling the ocean of the landsea mask.
    Ensures that region that may be cut-off from the open ocean due to very
    coarse resolution targets are filled using pre-defined seed points.

    Parameters
    ----------
    landseamask : :class:`iris.cube.Cube`
        Landsea mask.

    Returns
    -------
    ocean: :class:`iris.cube.Cube`
        Open ocean mask.

    """

    def get_seed_index(lat_lon):
        # Get seed point in index points
        (lat, lon) = lat_lon
        x, y = ants.utils.cube.horizontal_grid(landseamask, dim_coords=True)
        if (
            lon < x.points.min()
            or lon > x.points.max()
            or lat < y.points.min()
            or lat > y.points.max()
        ):
            # Ocean floodfill seed value is not contained within the domain
            # extent.
            return None
        xi = abs(x.points - lon).argmin()
        yi = abs(y.points - lat).argmin()

        xd = landseamask.coord_dims(x)[0]
        yd = landseamask.coord_dims(y)[0]
        seed = [0, 0]
        seed[yd], seed[xd] = yi, xi
        return seed

    # Define lat-lon ocean seed points for filling.
    seeds = [
        [35.1, 17.9],
        [43.3, 31.5],
        [83, 144],
        [59.49, 273.9],
        [14.26, 280.75],
        [25.59, 269.62],
        [52.79, 149.47],
        [39.95, 134.88],
        [-5.148, 112.34],
        [-6.149, 126.37],
        [0.655, 124.8],
        [39.23, 51.51],
    ]
    src_crs = ccrs.Geodetic()
    tgt_crs = landseamask.coord(axis="x").coord_system
    tgt_crs = tgt_crs.as_ants_crs().as_cartopy_crs()
    for seed in seeds:
        seed = tgt_crs.transform_point(seed[1], seed[0], src_crs)[::-1]
        seedi = get_seed_index(seed)
        if (
            seedi is None
            or landseamask.data[tuple(seedi)] == 1
            or landseamask.data[tuple(seedi)] == 2
        ):
            # Skip seed point as it is land (1), already filled (2) or isn't
            # contained within extent.
            continue
        ants.analysis.floodfill(landseamask.data, seedi, 2, extended_neighbourhood=True)
    ocean = landseamask.copy((landseamask.data == 2))
    return ants.utils.cube.defer_cube(ocean)


def make_consistent_with_mask(cube, ocean_mask):
    """
    Make cube consistent with the supplied ocean mask.

    Parameters
    ----------
    cube : :class:`iris.cube.Cube`
        Cube to made consistent with supplied ocean mask.
    ocean_mask : :class:`iris.cube.Cube`
        Target ocean mask.

    Returns
    -------
    : None
        In place operation

    """
    ocean = ocean_mask.data != 0
    cube.data[ocean] = 0


def get_total_gradient(orography):
    """
    Compute absolute gradients of supplied orography field using least
    squares of 2-D gradient components. Original blending code computed slightly
    differently. See get_total_gradient_original()
    """
    grad_x, grad_y = ants.analysis.calc_grad(orography)
    grad = (grad_x.data ** 2 + grad_y.data ** 2) ** 0.5

    return grad


def get_total_gradient_original(orography):
    """
    Method used to compute gradients of orography in code prior to
    the current version, retained for legacy purposes. Not currently
    used (ANTS 1.2)
    """
    radius = 6371229.0
    nrows = orography.shape[0]  # Latitudes
    max_row = 80.0  # Do not calculate dx poleward of this latitude (degrees)
    row_length = orography.shape[1]  # Longitudes
    dlambda = 2 * np.pi * radius / row_length
    dphi = np.pi * radius / nrows
    min_dx = np.cos(np.radians(max_row)) * dlambda
    dx = np.sin(((np.arange(0, nrows) + 0.5) / nrows) * np.pi) * dlambda
    dx = np.array([dx] * row_length).T
    dx[dx < min_dx] = min_dx
    dx[dx > dphi] = dphi
    grad_x, grad_y = np.gradient(orography.data)
    grad = (grad_x ** 2 + grad_y ** 2) ** 0.5
    grad /= dx
    return grad


def butterfilter(data, wavenum):
    """
    filter 2-D input data
    with filter order "order", wavenum=1/n (n samples per half-cycle)
    """
    order = 2
    Bconst, Aconst = signal.butter(order, wavenum, "high", output="ba")
    xdim = data.shape[1]
    ydim = data.shape[0]
    datafft = np.fft.fft2(data)
    freqx = 2 * np.pi * np.fft.fftfreq(xdim)
    freqy = 2 * np.pi * np.fft.fftfreq(ydim)  # noqa: E741
    freqxmesh, freqymesh = np.meshgrid(freqx, freqy)
    freq = (freqxmesh ** 2 + freqymesh ** 2) ** 0.5
    freqs, response = signal.freqz(Bconst, Aconst, worN=freq)
    datafft *= np.abs(response)
    datafilt = np.fft.ifft2(datafft)
    datalow = data - datafilt
    return np.real(datafilt), np.real(datalow)


def blended_orography(
    fine_orography,
    smooth_orography,
    max_slope,
    max_slope_from_smooth,
    max_smooth,
    dilation=False,
):
    smooth_grad = get_total_gradient(smooth_orography)

    # Set max_slope to that of smoothed orography, if desired
    if max_slope_from_smooth:
        max_slope = np.max(smooth_grad)

    fine_grad = get_total_gradient(fine_orography)
    isteep = np.where(fine_grad > max_slope)

    # Number of points for dilation and butterfilter
    # If global, use 120 km lengthscale
    if ants.utils.cube.is_global(fine_orography) and dilation:
        radius = fine_orography.coord(axis="y").coord_system.semi_major_axis
        dx = radius * np.diff(np.radians(fine_orography.coord(axis="y").points)).max()
        npts = np.max([int(120000.0 / dx), 2.0])
    else:
        npts = 5

    # Initialise loop variables
    smooth_count = 0
    print("Number of points to smooth =", isteep[0].size)
    # Keep looping until there are no blended slopes steeper than the "smooth" orog file
    while isteep[0].size > 0:

        # Calculate fine data stats
        print("Loop number = ", smooth_count)

        isteep = np.where(fine_grad > max_slope)

        targmatrix = np.zeros_like(smooth_grad)
        targmatrix[isteep] = 1

        # Perform dilation to expand region
        if dilation:
            targmatrix = ndimage.morphology.binary_dilation(
                targmatrix, iterations=int(npts * 2)
            ).astype(float)
        # low pass filter, limit max=1
        matrixhipass, targmatrix = butterfilter(targmatrix, 1.0 / npts)
        targmatrix[targmatrix > 1] = 1
        targmatrix[targmatrix < 0] = 0

        # blend the smooth and fine ancils
        blended = smooth_orography * targmatrix + fine_orography * (1 - targmatrix)

        fine_orography = blended.copy()

        fine_grad = get_total_gradient(fine_orography)

        isteep = np.where(fine_grad > max_slope)

        smooth_count += 1

        if smooth_count > max_smooth:
            print("We have smoothed over ", max_smooth, " times - break loop")
            break

    return fine_orography
