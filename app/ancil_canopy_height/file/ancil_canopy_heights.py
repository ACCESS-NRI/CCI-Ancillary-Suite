#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
"""
Canopy Heights
**************

Derive the canopy heights from the leaf area index by the following relation::

    Canopy height = height factor * LAI^(2/3)

Following this, override all tree PFTS using a trees dataset by the following
method:

- Regrid the tree source dataset to the target grid.
- Make this field consistent with the landsea mask of the LAI by performing an
  indexed based nearest neighbour search, constraining y to a distance of 50km
  (100km latitude band).
- Override tree fields we calculated from the above relation with those derived
  from this tree source.

The canopy heights returned represents the maximum value across the year
assigned to each month, due to the current basic handling of seasonal
variability.

"""
import iris

import ants
import ants.decomposition as decomp
import ants.io.save as save
import canopy_heights


def load_data(source_path, trees_filepath, canopy_height_factors_path):
    lai_stash_con = iris.AttributeConstraint(STASH='m01s00i217')
    #lai_cube = ants.io.load.load_cube(source_path, lai_stash_con)
    lai_cube = ants.io.load.load_cube(source_path, lai_stash_con)
    canopy_height_factors = canopy_heights.load(canopy_height_factors_path)
    trees_cube = ants.io.load.load_cube(trees_filepath)
    #trees_cube = ants.io.load.load_cube(trees_filepath)
    return lai_cube, trees_cube, canopy_height_factors


def duplicate_canopy_heights_over_12_months(canopy_heights_cube, lai):
    """
    Duplicate canopy heights across each month.

    Currently the canopy heights (which is a single time field) ends up
    stored in the same file as the periodic 12-month LAI fields (so dummy
    fields are necessary to do this).

    """
    new_canopy_heights = lai.copy()
    tdim = new_canopy_heights.coord_dims(new_canopy_heights.coord('time'))
    pdim = new_canopy_heights.coord_dims(
        new_canopy_heights.coord('pseudo_level'))
    ydim = new_canopy_heights.coord_dims(new_canopy_heights.coord(axis='y'))
    xdim = new_canopy_heights.coord_dims(new_canopy_heights.coord(axis='x'))
    if not (len(tdim) == len(pdim) == len(ydim) == len(xdim) == 1):
        msg = 'Expecting 1D dimension coordinates.'
        raise RuntimeError(msg)
    new_canopy_heights.data.transpose(
        tdim[0], pdim[0], ydim[0], xdim[0])[:] = (
            canopy_heights_cube.data[:])
    new_canopy_heights.metadata = canopy_heights_cube.metadata
    return new_canopy_heights


def main(source_filepath, trees_filepath, canopy_height_factors,
         tree_ids, output_filepath, netcdf_only):
    lai, trees_cube, canopy_height_factors = load_data(
        source_filepath, trees_filepath, canopy_height_factors)
    canopy_heights_cube = canopy_heights.calc_canopy_heights(
        lai, canopy_height_factors)

    # Apply the Simard trees dataset to our canopy heights.
    trees = decomp.decompose(ants.analysis.mean, trees_cube,
                             canopy_heights_cube[0])
    canopy_heights_cube = canopy_heights.calc_tree_pfts(
        trees, canopy_heights_cube, tree_ids)

    canopy_heights_cube = duplicate_canopy_heights_over_12_months(
        canopy_heights_cube, lai)
    
    save.netcdf(canopy_heights_cube, output_filepath)
    if not netcdf_only:
        save.ancil(canopy_heights_cube, output_filepath)
        
    return canopy_heights_cube


if __name__ == '__main__':
    parser = ants.AntsArgParser()
    parser.add_argument('--trees-dataset', type=str, required=True)
    parser.add_argument('--canopy-height-factors', type=str,
                        help='Canopy height factors.',
                        required=True)
    parser.add_argument('--tree-ids', type=str, default='1,101,102,103,2,201,202')
    args = parser.parse_args()
    tree_ids = [int(tree_id) for tree_id in args.tree_ids.split(',')]
    main(args.sources, args.trees_dataset, args.canopy_height_factors,
         tree_ids, args.output, args.netcdf_only)
