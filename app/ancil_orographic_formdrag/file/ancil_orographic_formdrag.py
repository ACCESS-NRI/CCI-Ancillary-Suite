#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
"""
Ororgraphic formdrag generation
*******************************
Generate the orographic formdrag ancillaries:

 - half of (peak to trough height of orography) (m01s00i018)
 - silhouette of orography per unit area (m01s00i017)

These fields are formed using the following steps:

- Regrid the provided silhouette and half height sources onto the target 
  grid and make the result consistent with the provided landseamask.
- Set all ocean points to zero, rather than leaving them set as missing data. 
  (Not doing this causes the UM to crash on the first time step.)
- Divide the data in the half height field by 2*sqrt(2) to match the variable 
  output by the CAP.

"""
import ants
import ants.decomposition as decomp
import ants.io.save as save
import iris
import numpy as np


def load_data(formdrag_filepath, target_filepath):
    half_height = ants.io.load.load_cube(formdrag_filepath, "m01s00i018")
    silhouette = ants.io.load.load_cube(formdrag_filepath, "m01s00i017")

    lsm_con = iris.AttributeConstraint(STASH="m01s00i030")
    target_lsm = ants.io.load.load_cube(target_filepath, constraint=lsm_con)
    ants.utils.cube.guess_horizontal_bounds(target_lsm)

    return silhouette, half_height, target_lsm


def apply_factor(half_height):
    """
    Divide the half height field by a factor of 2*sqrt(2).

    This division is performed in order to be consistent with the variable
    output by the CAP.

    Parameters
    ----------
    half_height : :class:`~iris.cube.Cube`
        Field representing half of the peak to trough height.

    Returns
    -------
    : None
        In-place operation

    """
    half_height.data = half_height.data / (2 * np.sqrt(2))


def set_sea_points(source, ocean):
    """
    Set the points denoted by the ocean array to zero, rather than leaving as
    missing data.

    Leaving the ocean points as missing data rather than setting them to zero
    causes the UM to crash on the first time step.

    Parameters
    ----------
    source : :class:`~iris.cube.Cube`
        Field representing either half of the peak to trough height or the
        silhouette.
    ocean : :class:`~numpy.ndarray`
        The array to denote the ocean points. Ocean points are set to True,
        all others to False.

    Returns
    -------
    : None
        In-place operation

    """
    source.data[ocean] = 0


def regrid_cube(source, target):
    scheme = ants.regrid.GeneralRegridScheme(horizontal_scheme="TwoStage")
    return source.regrid(target, scheme)


def main(sources, target, output, netcdf_only):
    silhouette, half_height, lbm = load_data(sources, target)

    # Regrid cubes
    silhouette_result = decomp.decompose(regrid_cube, silhouette, lbm)
    half_height_result = decomp.decompose(regrid_cube, half_height, lbm)

    ants.analysis.make_consistent_with_lsm(silhouette_result, lbm, True)
    ants.analysis.make_consistent_with_lsm(half_height_result, lbm, True)

    # Apply scaling factor to the half height cube to match the CAP output
    apply_factor(half_height_result)

    # Set ocean points to zero rather than missing data to avoid UM crash
    ocean = lbm.data == 0
    set_sea_points(silhouette_result, ocean)
    set_sea_points(half_height_result, ocean)

    save.netcdf([silhouette_result, half_height_result], output)
    if not netcdf_only:
        save.ancil([silhouette_result, half_height_result], output)


if __name__ == "__main__":
    parser = ants.AntsArgParser(target_lsm=True)
    args = parser.parse_args()
    main(
        args.sources,
        args.target_lsm,
        args.output,
        args.netcdf_only,
    )
