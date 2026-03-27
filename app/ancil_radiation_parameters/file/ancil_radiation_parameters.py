#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
"""
Orographic radiation parameters generation
******************************************
Calculate the unfiltered mean and gradients on the target grid. This processing
should be used in combination with a source orography that has not undergone the
Raymond filter preprocessing.

The following methodology is adopted:

- The mean orography is derived by conservatively regridding the source
  orographic height field onto the target grid.
- The grad x and grad y fields are calculated from the source orographic height
  field and conservatively regridded to the target grid. Any 'nan' points are
  set to 0.
- All fields are ensured to have values of 0 where the supplied landsea mask
  says there is open ocean. Open ocean is that region which is identified via a
  floodfill from the provided landsea mask.

This generates the following ancillaries:

 - unfiltered orography (m01s00i007)
 - slope_x (m01s00i005)
 - slope_y (m01s00i006)

"""
import functools

import ants
import ants.decomposition as decomp
import ants.io.save as save
import cartopy.crs as ccrs
import iris
import iris.analysis.calculus
import numpy as np
import orography_utils


def calc_slopes(source, ocean_mask):
    """
    Calculate the two gradient fields of the orography.

    Calculates the two gradient fields from the source orography. Regrids these
    two gradient fields to the target grid. Any points beyond the pole cut off
    (88.0 + 1e-10) are set to zero, as are any 'nan' points. Both fields are
    ensured to have values of 0 where the supplied landsea mask says there is
    open ocean. (Open ocean is that region which is identified via a floodfill
    from the provided landsea mask.) Also fixes the metadata (name, STASH and
    units) of the two fields.

    Parameters
    ----------
    source : :class:`iris.cube.Cube`
        Orography source dataset.
    ocean_mask : :class:`iris.cube.Cube`
        Target landsea mask.

    Returns
    -------
    : :class:`iris.cube.Cube`, :class:`iris.cube.Cube`
        slope_x (m01s00i005), slope_y (m01s00i006).

    """
    slope_x, slope_y = ants.analysis.calc_grad(source)

    lower_lats = source.coord(axis="y").bounds[:, 0]
    upper_lats = source.coord(axis="y").bounds[:, 1]

    pole_cut_off = 88.0 + 1e-10

    # Check for any points where the bounds exceed the pole cut off to avoid
    # missing points at low resolutions.
    slope_x.data[np.abs(lower_lats) > pole_cut_off, :] = 0.0
    slope_y.data[np.abs(lower_lats) > pole_cut_off, :] = 0.0

    slope_x.data[np.abs(upper_lats) > pole_cut_off, :] = 0.0
    slope_y.data[np.abs(upper_lats) > pole_cut_off, :] = 0.0

    regridder = ants.regrid.GeneralRegridScheme(horizontal_scheme="TwoStage").regridder(
        source, ocean_mask
    )

    slope_x = regridder(slope_x)
    slope_y = regridder(slope_y)

    slope_x.units = "m/m"
    slope_y.units = "m/m"
    slope_x.rename("Orographic Gradient X Component")
    slope_y.rename("Orographic Gradient Y Component")
    slope_x.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i005")
    slope_y.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i006")

    # Ensure consistency with the supplied ocean mask.
    for cube in [slope_x, slope_y]:
        orography_utils.make_consistent_with_mask(cube, ocean_mask)

    # Set nan values to 0.
    for cube in [slope_x, slope_y]:
        cube.data[np.isnan(cube.data)] = 0

    return (slope_x, slope_y)


def load_data(source, target):
    src_cube = ants.io.load.load_cube(source, "unfiltered_surface_altitude")
    tgt_cube = ants.io.load.load_cube(target, "land_binary_mask")
    ants.utils.cube.guess_horizontal_bounds(src_cube)
    ants.utils.cube.guess_horizontal_bounds(tgt_cube)
    return src_cube, tgt_cube


def main(sources, target, output, netcdf_only):
    orog_height, target = load_data(sources, target)

    mean = decomp.decompose(orography_utils.calc_mean, orog_height, target)

    ocean = orography_utils.gen_ocean_mask(target)

    orography_utils.make_consistent_with_mask(mean, ocean)

    mean.rename("unfiltered orography")
    mean.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i007")

    results_slopes = decomp.decompose(calc_slopes, orog_height, ocean)

    results = [mean] + results_slopes

    save.netcdf(results, output)
    if not netcdf_only:
        save.ancil(results, output)

    return results


if __name__ == "__main__":
    parser = ants.AntsArgParser(target_lsm=True)
    args = parser.parse_args()
    main(
        args.sources,
        args.target_lsm,
        args.output,
        args.netcdf_only,
    )
