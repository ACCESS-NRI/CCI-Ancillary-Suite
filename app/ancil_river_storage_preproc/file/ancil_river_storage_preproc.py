#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
Pre-processing river trip storage application
*********************************************

Corrects to the source data in order to allow its
interpretation.

The following changes are applied:

* Modify the source to represent a monthly mean, ensuring that a periodic
  timeseries can be recognised in the final ancillary.

"""
import warnings

import ants
import ants.config
import ants.fileformats.pp as pp
import ants.io.save as save
import ants.utils
import iris.fileformats.pp as ipp

STORAGE = {
    "stash_in": ipp.STASH(0, 26, 1),
    "stash_out": ipp.STASH(1, 0, 153),
    "long_name": "river_routing_storage",
}


def correct_metadata(filename):
    """Load the source data, correct the metadata and return Iris cubes."""
    fields = list(pp.load_ppfields(filename))
    fields = pp.field_filter(fields, STORAGE["stash_in"])

    def as_periodic_time_series(field):
        """
        Updates the field to make it look like a monthly mean.

        This then carries through to the final ancillary, setting it as
        periodic if source is a full 12 months (otherwise, identifies as time
        series).

        """
        field.lbtim.ia = 0
        field.lbtim.ib = 2
        field.lbtim.ic = 2
        field.lbdatd = 1
        # The data is representative of the 1950 - 2000 period.
        field.lbyr = 1950
        field.lbyrd = 1950
        if field.lbmond < 12:
            field.lbmond += 1
        else:
            field.lbmond = 1
            field.lbyrd += 1
        field.lbmin = 0
        field.lbmind = 0
        try:
            field.lbsec = 0
            field.lbsecd = 0
        except AttributeError:
            pass

    if not (fields[0].lbdatd == 0 and fields[0].lbtim.ib == 0):
        warnings.warn("Pre-processing metadata correction stage not required")
    else:
        for field in fields:
            as_periodic_time_series(field)

    storage_cube = pp.pp2cubes(fields)[0]
    storage_cube.long_name = STORAGE["long_name"]
    storage_cube.attributes["STASH"] = STORAGE["stash_out"]
    storage_cube.coord("time").attributes["representative_period"] = "1950 - 2000"

    return storage_cube


def main(source_filepath, output_filepath, use_new_saver, netcdf_only):
    cube = correct_metadata(source_filepath)

    if use_new_saver:
        save.netcdf(cube, output_filepath)
    else:
        ants.save(cube, output_filepath)

    return cube


def _get_parser():
    parser = ants.AntsArgParser()
    return parser


if __name__ == "__main__":
    args = _get_parser().parse_args()
    main(args.sources, args.output, args.use_new_saver, args.netcdf_only)
