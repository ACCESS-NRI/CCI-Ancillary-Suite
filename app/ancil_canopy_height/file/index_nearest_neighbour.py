# (C) Crown Copyright, Met Office. All rights reserved.
import cartopy.geodesic
import numpy as np


def ydist2index(cube, dist_limit):
    """
    Return the number of grid points in 'y' for the distance specified in 'y'.

    Assuming a regular grid, calculate the 'y' index distance corresponding to
    the specified 'y' distance in km.  Calculation is floating point based,
    where the index returned is rounding up to the nearest integer.

    Parameters
    ----------
    cube : :class:`~iris.cube.Cube`
        Source cube defining the grid no which to calculate the index-distance
        relation.
    dist_limit : tuple
        y distance limit in km.

    Returns
    -------
    : int
        y index distance.

    """
    def get_index(start_pnts, end_pnts, geodesic, dist):
        distance_per_lat_grid_pnt = geodesic.inverse(
            start_pnts, end_pnts).base[0, 0]

        return int(np.ceil((dist / distance_per_lat_grid_pnt)))

    dist_limit = np.array(dist_limit)
    # Convert 'km' units to 'm'
    dist_limit *= 1*(10**3)

    y = cube.coord(axis='y').points
    x = cube.coord(axis='x').points
    crs = cube.coord(axis='y').coord_system.as_ants_crs().as_cartopy_crs()
    print(f"crs is: {crs}")
    geodesic = cartopy.geodesic.Geodesic(crs.globe.semimajor_axis, 0)
    print(f"geodesic is: {geodesic}")

    # Gen y
    start_pnts = [x[0], y[0]]
    end_pnts = [x[0], y[1]]
    return get_index(start_pnts, end_pnts, geodesic, dist_limit)


def _neighbour_eval_2d(data, yy_xx, wraparound, loop_limit):
    shape = data.shape
    edge_length_y = edge_length_x = 1
    res = np.nan
    yy, xx = yy_xx
    seed_y, seed_x = yy, xx
    for lp in range(max(loop_limit)):
        if loop_limit[0] > lp:
            edge_length_y += 2
            seed_y = seed_y-1
        if loop_limit[1] > lp:
            edge_length_x += 2
            seed_x = seed_x-1
            if wraparound:
                seed_x %= shape[1]

        ypath = list(range(seed_y, seed_y + edge_length_y-1))  # left edge
        ypath += [ypath[-1] + 1] * (edge_length_x-1)  # top edge
        ypath += list(range(ypath[-1], ypath[0], -1))  # right edge
        ypath += [ypath[0]] * (edge_length_x-1)  # bottom edge
        ypath = np.array(ypath)

        xpath = [seed_x] * (edge_length_y-1)  # left edge
        xpath += list(range(xpath[-1], xpath[-1] + edge_length_x-1))  # top edge
        xpath += [xpath[-1] + 1] * (edge_length_y-1)  # right edge
        xpath += list(range(xpath[-1], xpath[0], -1))  # bottom edge
        xpath = np.array(xpath)

        if wraparound:
            xpath %= shape[1]
        # Exclude indices beyond array bounds.
        inside = ((xpath <= (shape[1]-1)) * (xpath >= 0) *
                  (ypath <= (shape[0]-1)) * (ypath >= 0))
        xpath = xpath[inside]
        ypath = ypath[inside]

        debug = False
        if debug:
            dbg_data = np.zeros(data.shape)
            dbg_data[yy, xx] = 3
            dbg_data[seed_y, seed_x] = 2
            dbg_data[ypath, xpath] = 1
            print('missing: 3, seed: 2, loop: 1')
            print(dbg_data)

        # Fill with mean value
        mean = data[ypath, xpath].mean()
        if not np.ma.is_masked(mean):
            res = mean
            break
    return res


def fill_nearest_index(data, target_mask, wraparound=False, loop_limit=None):
    """
    Fill nearest index data using index space.

    Iteratively expand our index based neighbourhood search until we find
    data to fill our missing data points.  Where more than one point is found
    within the same index distance, then the mean of these values is returned.

    Parameters
    ----------
    data : :class:`numpy.ma.ndarray`
        2D masked array to fill.
    target_mask : :class:`numpy.ndarray`
        2D bool array defining the target mask.  The difference between the
        provided source and the target mask is what is considered 'missing' to
        fill.
    wrapraround : bool, optional
        Search wraparound.
    loop_limit : tuple
        Constrain the iterative expansion of the index-based search
        neighbourhood in (y, x).  A value of -1 is used to denote unconstrained
        search in that given dimension (e.q. (6, -1) means to constrain 'y' to
        6 iterations (distance of 6 indices in y) while 'x' is unlimited
        (limited by only the shape of the provided data).

    Returns
    -------
    : None
        In-place operation.

    """
    # Missing values are those which are masked in data which are not masked in
    # the target.
    missing = data.mask * np.logical_not(target_mask)
    missing_indx = np.array(np.where(missing))
    replacement_value = np.ones(missing_indx.shape[1], dtype='float') * -1

    max_loop_limit = np.array(data.shape) - 1
    if loop_limit is None:
        loop_limit = max_loop_limit
    else:
        loop_limit = np.min([loop_limit, max_loop_limit], axis=0)
        max_loop_request = loop_limit == -1
        loop_limit[max_loop_request] = max_loop_limit[max_loop_request]
        loop_limit = np.max([loop_limit, [1, 1]], axis=0)

    for i in range(missing_indx.shape[1]):
        replacement_value[i] = _neighbour_eval_2d(data, missing_indx[:, i],
                                                  wraparound, loop_limit)
    data[missing_indx[0, :], missing_indx[1, :]] = replacement_value
    # Now inherrit the target mask.
    data.mask = target_mask.astype('bool', copy=False)
