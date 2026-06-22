#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
"""
Soil albedo generation
**********************
Generate the soil albedo (m01s00i216), which consists of:

- Regridding to the provided grid (lct).
- Making the albedo consistent with the land cover type fraction.
    - Setting soil albedo to 0.75 where the land cover type fraction says there
      is ice.
    - Utilise the spiral search to remedy differences between coastlines
      between the generated soil albedo and the land cover type fraction.

"""
import iris
import numpy as np

import ants
import ants.decomposition as decomp
import ants.io.save as save
import ants.utils

iris.FUTURE.save_split_attrs = True

def load_data(albedo_filepath, lct_filepath):
    albedo_cube = ants.io.load.load_cube(albedo_filepath)

    stash_con = iris.AttributeConstraint(STASH="m01s00i216")
    lct = ants.io.load.load_cube(lct_filepath, stash_con)
    ants.utils.cube.guess_horizontal_bounds(lct)
    return albedo_cube, lct


def regrid(source, target):
    return source.regrid(target, ants.regrid.GeneralRegridScheme(horizontal_scheme='TwoStage'))


def make_consistent_lct(soil_albedo, lct, ice_tile_id=9):
    """
    Make the soil albedo consistent with the provided land cover type fraction.

    Consistency involves:
    - Setting soil albedo where the lct identified as ice to 0.75
    - Filling through nearest neighbour remaining differences between soil
      albedo and lct coastlines.

    """

    def _fetch_lct(source, um_tile_id):
        pseudo_level = source.coord("pseudo_level")
        pseudo_level_points = pseudo_level.points.tolist()
        return pseudo_level_points.index(um_tile_id)

    # Set ice albedo value
    pseudo_level = lct.coord("pseudo_level")
    pdim = lct.coord_dims(pseudo_level)[0]
    ice_value = 0.75
    ice_level = _fetch_lct(lct, ice_tile_id)
    ice_soil_slices = [slice(None)] * lct.ndim
    ice_soil_slices[pdim] = ice_level
    ice = lct[tuple(ice_soil_slices)]
    soil_albedo.data[(ice.data >= 0.5) & ~np.ma.getmaskarray(ice.data)] = ice_value

    # We cannot ensure that soil albedo values are consistent with soil vs non
    # soil in the land cover type fraction as we don't know which soil albedo
    # values associate with what surface types...

    # Ensure ice is masked
    ants.utils.cube.fix_mask(ice)

    mask = ice.copy(ice.data.mask)
    filler = ants.analysis.FillMissingPoints(soil_albedo, target_mask=mask)
    filler(soil_albedo)
    return soil_albedo


def main(source_path, lct_filepath, ice_tile_id, output_filepath, netcdf_only):
    albedo, lct = load_data(source_path, lct_filepath)
    grid = lct[0]
    albedo = decomp.decompose(regrid, albedo, grid)
    make_consistent_lct(albedo, lct, ice_tile_id)

    save.netcdf(albedo, output_filepath)
    if not netcdf_only:
        save.ancil(albedo, output_filepath)

    return albedo


if __name__ == "__main__":
    parser = ants.AntsArgParser()
    parser.add_argument(
        "--lct-ancillary",
        type=str,
        help=("Filepath to the land cover type fraction " "ancillary."),
        required=True,
    )
    parser.add_argument(
        "--ice-tile-id",
        type=int,
        help="ID for the ice tile",
        default=9
    )
    args = parser.parse_args()
    main(args.sources, args.lct_ancillary, args.ice_tile_id, args.output, args.netcdf_only)
