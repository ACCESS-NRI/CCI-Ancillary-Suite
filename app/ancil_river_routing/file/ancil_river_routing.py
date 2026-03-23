#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
River routing application
*************************

Derive the river routing sequence and direction ancillaries.

The following steps are taken:

* Landcover fraction field (LCF) is regridded to the river routing field
  using an area weighted regrid.
* Regridded LCF is used to identify points which are entirely ocean and
  those coastal points which do not have any river flow going into them.
* At points that are entirely ocean and coastal points that have no river
  flow going into them, the river direction and sequence is set to missing
  data.
* At coastal points that do have flow going into them, the river direction
  is set to 9, to signify a pour point into the ocean.
* order_nemo_rivers: Sort rivers into an order, giving each river outflow
  a unique ID number. Puts these river numbers into UM and NEMO
  river_number ancillary files.
* On output of ancillary, the UM grid definition is inherited from the
  input LCF, as there is a detachment between the river routing field and
  the UM grid definition.

In the resulting river routing ancillary file the directions define where
the flow goes from point x:

 | 8  1  2
 | 7  x  3
 | 6  5  4

  9 = outflow to ocean
  10 = inland basin flow

e.g. 3 means the river flows eastwards

For GC5 (and later) ancillary files: "land_threshold_outflow=0.5" should be
provided to enable the outflow points to fall on the coast.
"make_nemo_rivers" should be set to True to enable the generation of UM and
NEMO river_number ancillary files.

For GC4 (and earlier) ancillary files: "land_threshold_outflow=0.0" and
"make_nemo_rivers=False" should replicate old river coupling behviour.

"""
import ants
import ants.config
import ants.decomposition as decomp
import ants.io.save as save
import ants.utils
import iris
import numpy as np
from proc_ants import order_nemo_rivers

POUR_POINT_INDICATOR = 9


def flow_in(direction_cube):
    """
    Determine which grid points have river flow going into them

    Parameters
    ----------
    direction_cube : :class:`iris.cube.Cube`
       River direction as output from ancil_river_routing.river_routing

    Returns
    -------
    flow_in_array : :class:`numpy.ndarray`
       Array of bool where True = Flow is going into this point and
       False = Flow is not going into this point

    """

    # Roll the cubes to get values in all directions
    west = np.roll(direction_cube.data, 1, axis=1)
    north = np.roll(direction_cube.data, 1, axis=0)
    east = np.roll(direction_cube.data, -1, axis=1)
    south = np.roll(direction_cube.data, -1, axis=0)
    northwest = np.roll(west, 1, axis=0)
    northeast = np.roll(east, 1, axis=0)
    southwest = np.roll(west, -1, axis=0)
    southeast = np.roll(east, -1, axis=0)

    # Flow is direction it is going in (e.g. 3 = eastwards)
    # 8  1  2
    # 7  9  3
    # 6  5  4
    # 9 = outflow to ocean
    # 10 = inland basin flow

    # See if any of the neighbours flow into this grid box
    flow_in_array = (
        (west == 3)
        | (north == 5)
        | (east == 7)
        | (south == 1)
        | (northwest == 4)
        | (northeast == 6)
        | (southwest == 2)
        | (southeast == 8)
    )
    return flow_in_array


def fix_nans(lcf_cube):
    """
    Sometimes decomp.decompose returns NANs along the northern and southern rows.
    We know that there is land on the southern row (Antarctica) and sea on the
    northern row (Arctic Ocean) so we should fix this here.

    Parameters
    ----------
    lcf_cube : :class:`iris.cube.Cube`
       Land cover fraction on the TRIP grid

    """

    # Make northern NAN points sea
    where_nan = np.where(np.isnan(lcf_cube.data[0, :]))
    if len(where_nan[0]) > 0:
        lcf_cube.data[0, where_nan] = 0.0

    # Make southern NAN points land
    where_nan = np.where(np.isnan(lcf_cube.data[-1, :]))
    if len(where_nan[0]) > 0:
        lcf_cube.data[-1, where_nan] = 1.0


def load_data(source, land_cover_fraction):
    sequence_cube, direction_cube = ants.io.load.load_cubes(
        source, ["river_routing_sequence", "river_routing_direction"]
    )
    lcf_cube = ants.io.load.load_cube(land_cover_fraction, "land_area_fraction")

    return sequence_cube, direction_cube, lcf_cube


def river_routing(
    sequence_cube, direction_cube, land_cover_fraction_cube, land_threshold_outflow=0.0
):
    # Regrid the land cover fraction to the river trip routing source grid.
    lcf_cube = decomp.decompose(
        ants.analysis.mean, land_cover_fraction_cube, direction_cube
    )

    fix_nans(lcf_cube)

    # Find out what is mostly sea and what is all sea
    mostly_sea = lcf_cube.data <= land_threshold_outflow
    all_sea = lcf_cube.data == 0.0

    # Set pour point where its next to the coast.
    direction_cube.data[mostly_sea] = POUR_POINT_INDICATOR

    # Mask sea locations where there is no flow going into the sea
    no_flow_in = ~flow_in(direction_cube)
    mask = all_sea & no_flow_in
    direction_cube.data = np.ma.array(direction_cube.data, mask=mask, copy=False)
    sequence_cube.data = np.ma.array(sequence_cube.data, mask=mask, copy=False)

    return sequence_cube, direction_cube


def main(
    source_filepath,
    land_cover_fraction_filepath,
    output_filepath,
    make_nemo_rivers,
    ocean_runoff_file,
    um_river_number_file,
    ocean_river_number_file,
    orca_dom_file,
    land_threshold_outflow,
    netcdf_only,
):
    sequence_cube, direction_cube, lcf_cube = load_data(
        source_filepath, land_cover_fraction_filepath
    )
    sequence_cube, direction_cube = river_routing(
        sequence_cube,
        direction_cube,
        lcf_cube,
        land_threshold_outflow=land_threshold_outflow,
    )

    if make_nemo_rivers:
        # Note order is (direction, sequence, ...) while most arguments in
        # this file are (sequence, direction, ...)
        order_nemo_rivers.main(
            direction_cube,
            sequence_cube,
            ocean_runoff_file=ocean_runoff_file,
            orca_dom_file=orca_dom_file,
            um_river_number_file=um_river_number_file,
            ocean_river_number_file=ocean_river_number_file,
        )

    # Note order is (direction, sequence, ...) while most arguments in
    # this file are (sequence, direction, ...)
    cubes = iris.cube.CubeList([direction_cube, sequence_cube])
    
    save.netcdf(cubes, output_filepath)
    if not netcdf_only:
        save.ancil(cubes, output_filepath)

    return sequence_cube, direction_cube


def _get_parser():
    parser = ants.AntsArgParser()
    parser.add_argument(
        "--land-cover-fraction", type=str, help="Land cover fraction.", required=True
    )
    parser.add_argument(
        "--make_nemo_rivers",
        help="Make NEMO river number ancillary file",
        action="store_true",
        required=False,
        default=False,
    )
    parser.add_argument(
        "--ocean_runoff_file",
        type=str,
        help="Ocean runoff used in ocean only runs.",
        required=False,
    )
    parser.add_argument(
        "--um_river_number_file",
        type=str,
        help="Output file for UM river number.",
        required=False,
    )
    parser.add_argument(
        "--ocean_river_number_file",
        type=str,
        help="Output file for NEMO river number.",
        required=False,
    )
    parser.add_argument(
        "--orca_dom_file",
        type=str,
        help="NEMO domain file (contains land-sea mask for NEMO model on ORCA grid)",
        required=False,
    )
    parser.add_argument(
        "--land_threshold_outflow",
        help="Land fraction threshold below which rivers outflow into the sea. "
        + "Default is 0 (i.e. needs to be all sea before the river arrives at "
        + "an outflow point).",
        type=float,
        required=False,
        default=0.0,
    )
    return parser


if __name__ == "__main__":
    args = _get_parser().parse_args()

    main(
        source_filepath=args.sources,
        land_cover_fraction_filepath=args.land_cover_fraction,
        output_filepath=args.output,
        make_nemo_rivers=args.make_nemo_rivers,
        ocean_runoff_file=args.ocean_runoff_file,
        um_river_number_file=args.um_river_number_file,
        ocean_river_number_file=args.ocean_river_number_file,
        orca_dom_file=args.orca_dom_file,
        land_threshold_outflow=args.land_threshold_outflow,
        netcdf_only=args.netcdf_only,
    )
