#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
River trip storage application
******************************

Derive the river routing ancillaries by utilising the input storage and river
direction data.

The following steps are required in the derivation:

* The river storage source is first regridded via area weighted mean
  to the same grid as the direction ancillary.
* Mask status of elements of the direction cube are applied to the storage
  cube.

"""
import ants
import ants.io.save as save
import ants.utils
import iris
import numpy as np


def load_data(source_path):
    direction_constraint = iris.AttributeConstraint(STASH="m01s00i152")
    direction_cube = ants.load_cube(source_path, direction_constraint)
    storage_constraint = iris.AttributeConstraint(STASH="m01s00i153")
    storage_source = ants.load_cube(source_path, storage_constraint)
    return direction_cube, storage_source


def river_storage(direction_cube, storage_source):
    storage_cube = ants.analysis.mean(storage_source, direction_cube)
    storage_cube.data = np.ma.array(storage_cube.data, copy=False)
    # Ensure direction_cube has a mask
    ants.utils.cube.fix_mask(direction_cube)
    storage_cube.data[:, direction_cube.data.mask] = np.ma.masked
    return storage_cube


def main(source_path, output_filepath, use_new_saver, netcdf_only):
    direction_cube, storage_cube = load_data(source_path)
    storage_cube = river_storage(direction_cube, storage_cube)

    if use_new_saver:
        if not netcdf_only:
            save.ancil(storage_cube, output_filepath)
        save.netcdf(storage_cube, output_filepath)
    else:
        ants.save(storage_cube, output_filepath)

    return storage_cube


def _get_parser():
    parser = ants.AntsArgParser()
    return parser


if __name__ == "__main__":
    args = _get_parser().parse_args()
    main(args.sources, args.output, args.use_new_saver, args.netcdf_only)
