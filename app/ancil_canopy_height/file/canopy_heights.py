# (C) Crown Copyright, Met Office. All rights reserved.
import warnings

import iris
import numpy as np

import ants
import ants.fileformats.json
import ants.utils
import index_nearest_neighbour


def calc_tree_pfts(trees, canopy_heights, tree_ids = [1, 101, 102, 103, 2, 201, 202]):
    """
    Inject tree PFTs from the trees dataset.

    - Regrid the tree source dataset to the target grid.
    - Override the existing tree fields with those derived from the tree
      source.

    Parameters
    ----------
    trees : :class:`iris.cube.Cube`
        Trees dataset.
    canopy_heights : :class:`iris.cube.Cube`
        Canopy heights as derived from the LAI.

    Returns
    -------
    None (in-place operation)

    """
    # Ensure tree PFT masks are consistent with the lsm
    ants.utils.cube.fix_mask(canopy_heights)
    target_mask = canopy_heights[0].data.mask

    loop_lim_y = index_nearest_neighbour.ydist2index(trees, 500)
    if np.ma.is_masked(trees.data):
        index_nearest_neighbour.fill_nearest_index(
            trees.data,
            target_mask,
            loop_limit=(loop_lim_y, -1),
            wraparound=trees.coord(axis="x").circular,
        )

    # The following masking fix seems to be made redundant by #1866.
    # #1917 has been created  to resolve this.
    if np.ma.isMaskedArray(trees.data) and not np.ma.isMaskedArray(canopy_heights.data):
        canopy_heights.data = np.ma.array(canopy_heights.data)

    # Override canopy heights tree PFTs.
    for tree_id in tree_ids:
        try:
            tree_slices = ants.analysis.cover_mapping.fetch_lct_slices(
                canopy_heights, tree_id
            )
            canopy_heights.data[tree_slices] = trees.data
        except ValueError as err:
            message = "{} is not in list".format(tree_id)
            if message in err.args[0]:
                # Ignore case where the PFT is missing.
                pass
            else:
                raise err
    return canopy_heights


def calc_canopy_heights(lai, canopy_height_factors):
    """
    Calculate the canopy heights.

    Calculating the canopy heights is achieved by the following equation:
    ::
        canopy_height[PFT] = height_factor[PFT] * LAI[PFT]^(2/3)
    by the height factors, one for each of the JULES classes.

    The canopy heights for each month represents the maximum across the year,
    due to the current basic handling of seasonal variability.

    Parameters
    ----------
    lai : :class:`iris.cube.Cube`
        Leaf area index ancillary.
    canopy_height_factors : :class:`CanopyHeightFactor`
        Canopy height factors for the given JULES classes.

    Returns
    -------
    : :class:`iris.cube.Cube`
        Canopy heights.

    """
    # Extract functional types in the correct order.
    p_coord = lai.coord("pseudo_level")
    pdim = lai.coord_dims(p_coord)
    if len(pdim) != 1:
        raise ValueError(
            "Expecting a 1-dimensional mapping of pseudo-levels "
            "within the LAI, got: {}".format(len(pdim))
        )
    pdim = pdim[0]

    # Ensure that the canopy height factors are ordered as they are in the
    # LAI source.
    canopy_height_factors.reorder_height_factors(p_coord.points)

    canopy_heights = lai.copy()
    data = canopy_heights.data
    data = data ** (2.0 / 3.0)
    canopy_heights.data = data

    transpose_index = list(range(canopy_heights.ndim))
    transpose_index.pop(transpose_index.index(pdim))
    transpose_index += [pdim]
    data = canopy_heights.data.transpose(transpose_index)
    data *= canopy_height_factors.height_factors

    # Tidy-up the resulting cubes metadata.
    canopy_heights.rename("canopy_height")
    canopy_heights.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi(
        "m01s00i218"
    )
    canopy_heights.units = "m"

    # Return the maximum Canopy heights over the year assigned to each of the
    # 12 months.
    tdim = canopy_heights.coord_dims(canopy_heights.coord("time"))
    if len(tdim) != 1:
        raise ValueError(
            "Expecting a 1-dimensional mapping of time "
            "within the canopy heights, got: {}".format(len(tdim))
        )
    tdim = tdim[0]
    slices = [slice(None)] * canopy_heights.ndim
    slices[tdim] = 0
    max_canopy_heights = canopy_heights[tuple(slices)].copy(
        np.amax(canopy_heights.data, axis=tdim)
    )
    max_canopy_heights.remove_coord("time")
    return max_canopy_heights


class CanopyHeightFactor(object):
    """Canopy height factor container."""

    def __init__(self, height_factors, jules_classes):
        self._height_factors = np.asarray(height_factors, dtype="float")
        self._jules_classes = np.asarray(jules_classes, dtype="int")
        self._valid_jules_classes()

    def _valid_jules_classes(self, jules_classes=None):
        """Check provided JULES_classes are valid."""
        if jules_classes is not None:
            jules_classes = np.asarray(jules_classes)
            if set(jules_classes) != set(self._jules_classes):
                msg = (
                    "The provided JULES Classes does not agree with those "
                    "mapped with the canopy height factors."
                )
                raise ValueError(msg)
        height_factors = self._height_factors
        jules_classes = (
            jules_classes if jules_classes is not None else self._jules_classes
        )

        if height_factors.ndim != 1 or jules_classes.ndim != 1:
            raise ValueError("Expecting 1D JULES classes and height factors.")
        if len(set(jules_classes)) != jules_classes.size:
            raise ValueError("Duplicate JULES classes specified")
        if height_factors.size != jules_classes.size:
            raise ValueError("Expecting 1 canopy height factor for each JULES " "class")

    @property
    def height_factors(self):
        """np.ndarray : Canopy height factors."""
        return self._height_factors

    @property
    def jules_classes(self):
        """np.ndarray : Jules classes."""
        return self._jules_classes

    def reorder_height_factors(self, jules_classes):
        """
        Re-order the height factors based on the JULES classes provided.

        Parameters
        ----------
        jules_classes : iterable of ints

        """
        jules_classes = np.asarray(jules_classes, dtype="int")
        self._valid_jules_classes(jules_classes)

        order = [
            np.argmax(self._jules_classes == jules_class)
            for jules_class in jules_classes
        ]
        self._height_factors = self._height_factors[order]
        self._jules_classes = jules_classes


def load(filename):
    """
    Load the canopy heights from a json file on disk.

    Parameters
    ----------
    filename : str
        Filepath to the json canopy height factors.  json file should contain
        both "canopy_heights" and "jules_classes" keys.

    Returns
    -------
    : CanopyHeightFactor
        Container object for the canopy height factors.

    ..warning:
        The form in which this information should take is still not known.  In
        practice this means the functionality or presence of this function may
        change.

    .. note:: Deprecated in 0.2dev:
        Function is to be replaced.  Unknown replacement at this time.

    """
    msg = "canopy_heights.load is experimental and may be deprecated."
    warnings.warn(msg, FutureWarning)

    required = ["canopy_height_factors", "jules_classes"]
    res = ants.fileformats.json.load(filename, keys=required, dtypes=["float", "int"])
    return CanopyHeightFactor(res["canopy_height_factors"], res["jules_classes"])
