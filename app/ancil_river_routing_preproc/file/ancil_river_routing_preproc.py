#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
Pre-processing river trip routing application
*********************************************

Apply selected modifications to the river routing dataset.

The following changes are applied:

* Modify the source time points so that they do not exceed 30 days in a
  month which would otherwise be invalid for a 360-day model calendar.

"""
import warnings

import ants
import ants.config
import ants.fileformats.pp as pp
import ants.io.save as save
import ants.utils
import iris
import iris.fileformats.pp as ipp

RIVER_SEQUENCE = {"stash": ipp.STASH(1, 0, 151), "long_name": "river_routing_sequence"}


RIVER_DIRECTION = {
    "stash": ipp.STASH(1, 0, 152),
    "long_name": "river_routing_direction",
}


def load_data(filename):
    fields = list(pp.load_ppfields(filename))

    direction = pp.field_filter_strict(fields, RIVER_DIRECTION["stash"])
    sequence = pp.field_filter_strict(fields, RIVER_SEQUENCE["stash"])
    return direction, sequence


def correct_metadata(direction, sequence):
    """Correct pp field metadata and return Iris cubes."""
    if not (direction.lbdatd == 31 and direction.lbtim == 2):
        warnings.warn("Pre-processing metadata correction stage not required")
    else:
        sequence.lbdatd = 30
        direction.lbdatd = 30

    sequence_cube = pp.pp2cubes(sequence)[0]
    sequence_cube.long_name = RIVER_SEQUENCE["long_name"]
    direction_cube = pp.pp2cubes(direction)[0]
    direction_cube.long_name = RIVER_DIRECTION["long_name"]

    return iris.cube.CubeList([sequence_cube, direction_cube])


def main(source_filepath, output_filepath, use_new_saver, netcdf_only):
    direction_field, sequence_field = load_data(source_filepath)
    cubes = correct_metadata(direction_field, sequence_field)

    if use_new_saver:
        save.netcdf(cubes, output_filepath)
    else:
        ants.save(cubes, output_filepath)

    return cubes


def _get_parser():
    parser = ants.AntsArgParser()
    return parser


if __name__ == "__main__":
    args = _get_parser().parse_args()
    main(args.sources, args.output, args.use_new_saver, args.netcdf_only)
