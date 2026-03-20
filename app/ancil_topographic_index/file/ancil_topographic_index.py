#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
Topographic index application
*****************************

Derive the mean and stdev topographic ancillaries and ensuring it's consistent
with the provided land cover type fraction ancillary.

The following steps are taken:

* Area weighted regrid the source to the grid specified by the
  land-sea-mask, see :func:`ants.analysis.mean`.
* Derive the standard deviation field by using the mean field just derived,
  see :func:`ants.analysis.stdev`.
* Update the mean and standard deviation field data so that they are consistent
  with the provided land cover type fraction fields (landsea mask and ice
  fraction fields).

  - Fill all non-ice locations with non-zero topographic index values using
    :func:`ants.analysis.FillMissingPoints`.
  - Set topographic index values to 0 for locations where there is ice.

.. note::
  GA7 generated ancillaries used a different algorithm for making the mean
  and standard deviation fields consistent with the land sea mask.  This
  was based on a Moore neighbourhood fill.

"""
import ants
import ants.config
import ants.decomposition as decomp
import ants.io.save as save
import iris
import numpy as np


def mean_stdev(source, target):
    mean_cube = ants.analysis.mean(source, target)
    stdev_cube = ants.analysis.standard_deviation(source, mean_cube)
    return mean_cube, stdev_cube


def load_data(source, lct):
    src_cube = ants.io.load.load_cube(source)
    lct_stash_con = iris.AttributeConstraint(STASH="m01s00i216")
    lct_cube = ants.io.load.load_cube(lct, lct_stash_con)
    ants.utils.cube.guess_horizontal_bounds([src_cube, lct_cube])
    return src_cube, lct_cube


def topographic_index(src_cube, lct_cube, ice_id):
    ice_level = ants.analysis.cover_mapping.fetch_lct_slices(lct_cube, ice_id)
    ice_level_cube = lct_cube[ice_level]
    target_mask = ice_level_cube.copy(np.ma.getmaskarray(ice_level_cube.data))

    result = decomp.decompose(mean_stdev, src_cube, ice_level_cube)
    # For now we write to separate files
    mean_cube = result.extract_cubes("mean topographic index")[0]
    mean_cube.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i274")

    # Don't fill locations which are supposed to be ice.  Topographic index
    # ice locations are masked.
    if not np.ma.isMaskedArray(ice_level_cube.data):
        ice_level_cube.data = np.ma.array(ice_level_cube.data, mask=False)

    ice = (ice_level_cube.data.data > 0) & (~np.ma.getmaskarray(ice_level_cube.data))
    target_mask_noice = mean_cube.copy(ice + target_mask.data)
    filler = ants.analysis.FillMissingPoints(mean_cube, target_mask=target_mask_noice)
    filler(mean_cube)
    stdev_cube = result.extract_cubes("standard deviation topographic index")[0]
    stdev_cube.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i275")
    filler(stdev_cube)

    # Lastly, ensure that ICE points are set to 0 rather than masked.
    mean_cube.data[ice] = 0
    stdev_cube.data[ice] = 0

    cubes = iris.cube.CubeList([mean_cube, stdev_cube])
    return cubes


def main(source_path, lct_path, output_filepath, ice_tile_id, netcdf_only):
    """
    Topographic index application top level call function.

    This top level function call loads the data, calls the underlying
    operations to derive the ancillaries and saves them to the requested
    location.

    Parameters
    ----------
    source_path : str
        Source data file path.
    lct_path : str
        File path to the land cover type fraction ancillary.
    output_filepath : str
        Output file path for the topographic index ancillary.

    Returns
    -------
    : tuple(:class:`~iris.cube.Cube`, :class:`~iris.cube.Cube`)
        Mean and standard deviation topographic index cubes.

    """
    source, lct = load_data(source_path, lct_path)
    cubes = topographic_index(source, lct, ice_tile_id)

    save.netcdf(cubes, output_filepath)
    if netcdf_only:
        save.ancil(cubes, output_filepath)

    return cubes


def _get_parser():
    parser = ants.AntsArgParser()
    lct_help = (
        "Land cover type fraction ancillary, used for ensuring "
        "consistency of the topographic index with the ice field and "
        "the landsea mask."
    )
    parser.add_argument("--lct-ancillary", type=str, help=lct_help, required=True)
    parser.add_argument("--ice-tile-id", type=int, help="ID of the ice tile", default=9)
    return parser


if __name__ == "__main__":
    args = _get_parser().parse_args()
    main(args.sources, args.lct_ancillary, args.output, args.ice_tile_id, args.netcdf_only)
