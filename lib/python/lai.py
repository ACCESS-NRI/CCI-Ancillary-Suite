# (C) Crown Copyright, Met Office. All rights reserved.
#
import warnings

import ants
import iris
import numpy as np


MIN_LAI = 1e-10


def fill_missing_data(lct, lai):
    """
    Populate the LAI missing data with appropriate values.

    In-place operation to fill data which is missing accross some months by
    setting to 1e-10 Set unmasked LAI < 1e-10 to 1e-10

    Parameters
    ----------
    lct : :class:`~iris.cube.Cube`
        Land cover type fraction, used to define the landsea mask.
        (PFT, Y, X).
    lai : :class:`~iris.cube.Cube`
        The LAI we want to process for missing data (TIME, Y, X).

    """
    dates = lai.coord('time').units.num2date(
        lai.coord('time').points)
    years = set([date.year for date in dates])
    months = sorted([date.month for date in dates])
    if len(years) != 1 or range(1, 13) != months:
        msg = ('Filling of LAI data expects a 12 month climatology for LAI '
               'total.\n{}')
        ValueError(msg.format(lai.coord('time')))

    # Apply the landsea mask
    pdim = lct.coord_dims(lct.coord('pseudo_level'))[0]
    lai.data[..., np.ma.getmaskarray(lct[pdim].data)] = np.ma.masked

    # Missing values across SOME months (i.e. values are missing across some
    # rather than all months).  Set these points to min_lai.
    tdim = lai.coord_dims(lai.coord('time'))[0]
    nn_missing = lai.data.mask.sum(axis=tdim)
    nn_missing = (nn_missing > 0) * (nn_missing < lai.shape[tdim])
    lai.data[lai.data.mask * nn_missing] = MIN_LAI

    # Those missing accross all months are handled when making the LAI
    # consistent with the coastlines (i.e. through nearest neighbour)

    # Set unmasked points with < minimum lai to the minimum LAI.
    lai.data[
        (lai.data.data <= MIN_LAI) * ~np.ma.getmaskarray(lai.data)] = MIN_LAI


def _calc_lai(lct, lai_total, relative_weights):
    """
    Calculate LAI by solving a set of linear equations through substitution.

    Where the total vegetation fraction (for the functional types) equals less
    than 10% of a grid box but non-zero, each LAI PFT is masked (to be filled
    using nearest neighbour when ensuring consistency with the lct coastlines).

    Parameters
    ----------
    lct : :class:`numpy.ndarray`
        The vegetation fraction array with order (pseudo-level, y, x).
    lai_total : numpy.ndarray
        The total lai array with order (t, y, x).
    relative_weights : :class:`numpy.ndarray`
        An array representing the relative weights between LAI classes.

    Returns
    -------
    : numpy.ndarray
        LAI result for each functional type (t, pseudo-level, y, x).

    """
    # Solve the set of linear equations using forward substitution.
    t1 = (relative_weights[:, None, None] * lct).sum(axis=0)
    t2 = relative_weights[:, None, None, None] * lai_total
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = t2 / t1
    if result.ndim != 4:
        msg = 'Unexpeced dimensionality: {}D, expecting {}D'
        raise ValueError(msg.format(result.ndim, 4))

    timedim = 1
    pldim = 0
    # Make our result (time, pseudo-level, y, x)
    result = result.transpose(timedim, pldim, 2, 3)
    # Set lai to zero where there is zero vegetation fraction.
    result[:, lct == 0] = 0
    # Set lai to zero where there is no total LAI.
    result.transpose(1, 0, 2, 3)[:, lai_total <= 0] = 0

    # Define a minimum veg fraction threshold (this ensures that the LAI
    # resolved doesn't give unreasonably high LAI PFTs in the case of a high
    # total LAI).  This threshold value was chosen by Martin Best.
    # We mask these elements which is then filled by nearest neighbour.
    min_veg_threshold = 0.2
    total_veg = lct.sum(axis=0)
    if not np.ma.is_masked(result):
        result = np.ma.array(result)
    result[..., (total_veg < min_veg_threshold) * (total_veg > 0)] = (
        np.ma.masked)

    # Ensure we retain the mask from the total LAI field (we may have been
    # overridden by the previous steps).
    if np.ma.is_masked(lai_total):
        result.mask.transpose(1, 0, 2, 3)[:, lai_total.mask] = True
    return result


def calc_lai(relative_weights, lai_total, lct):
    """
    Calculate the LAI for the vegetation fraction functional types.

    Parameters
    ----------
    relative_weights : dict
        Dictionary containing both weights ('relative_weights') and mappings
        ('jules_classes').
    lct : :class:`~iris.cube.Cube`
        Land cover type fraction ancillary containing the functional type
        classes.
    lai_total : :class:`~iris.cube.Cube`
        Total leaf area index data.

    Returns
    -------
    : :class:`~iris.cube.Cube`
        Lai corresponding to each functional type.

    """
    # Extract the horizontal grid definition from the land cover type fraction.
    slices = [0] * lct.ndim
    x, y = ants.utils.cube.horizontal_grid(lct)
    slices[lct.coord_dims(x)[0]] = slice(None)
    slices[lct.coord_dims(y)[0]] = slice(None)
    grid = lct[tuple(slices)]

    # Extract functional types.
    p_coord = lct.coord('pseudo_level')
    p_points = p_coord.points.tolist()
    slices = [slice(None)] * lct.ndim
    slices[lct.coord_dims(p_coord)[0]] = np.array(
        [p_points.index(val) for val in relative_weights['jules_classes']])

    slices = ants.analysis.cover_mapping.fetch_lct_slices(
        lct, relative_weights['jules_classes'])
    functional_types = lct[slices]

    # Regrid the LAI onto the lct grid.
    lai_tar = ants.analysis.mean(lai_total, grid)
    if not np.ma.isMaskedArray(lai_tar.data):                                  
        lai_tar.data = np.ma.array(lai_tar.data)
    fill_missing_data(functional_types, lai_tar)

    # Calculate the LAI for the functional types.
    results = _calc_lai(functional_types.data, lai_tar.data,
                        relative_weights['relative_weights'])

    # Tidy-up the resulting cubes metadata.
    lai_res = iris.cube.Cube(results)
    lai_res.standard_name = 'leaf_area_index'
    lai_res.units = 1
    lai_res.attributes['STASH'] = iris.fileformats.pp.STASH.from_msi(
        'm01s00i217')
    lai_res.add_dim_coord(y, 2)
    lai_res.add_dim_coord(x, 3)
    lai_res.add_aux_coord(functional_types.coord('pseudo_level'), 1)
    lai_res.add_dim_coord(lai_tar.coord('time'), 0)

    # Set unmasked points with < minimum lai to the minimum LAI.
    lai_res.data[(lai_res.data.data <= MIN_LAI) * ~lai_res.data.mask] = MIN_LAI

    return lai_res
