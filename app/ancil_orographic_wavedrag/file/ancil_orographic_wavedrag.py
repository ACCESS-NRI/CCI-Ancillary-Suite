#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
"""
Orographic wavedrag generation
******************************
Calculate the mean, standard deviation, sub-grid gradient correlation terms
and LM97 parameters on the target grid.

The following methodology is adopted:

- Calculate the filtered source orographic height (to remove scales < 6 km,
  computed in ancil_orographic_wavedrag_preproc.py) with mean gradient removed (SMG).
    - Linearly regrid the filtered mean above onto the source grid and remove
      this from the filtered source orography (surface_altitude_filtered, provided
      via the --filtered-source command line argument).
- Calculate the sub grid scale orography (SSO) fields:
    - Calculate grad x and grad y of the SMG and ensure that these are
      co-located.
        - Once calculated, grad x and grad y are staggered.  These are shifted
          by simple linear interpolation to co-locate them i.e.
          `p3 = (p1 + p2) / 2`.
    - Calculate the sub-grid correlation terms (grad_xx, grad_yy and grad_xy)
      and regrid them to the target grid.
    - Calculate the standard deviation field using the SMG and regridded SMG.
    - All fields are ensured to have values of 0 where the supplied landsea
      mask says there is open ocean.  Open ocean is that region which is
      identified via a floodfill from the provided landsea mask.
    - Set 'nan' to 0 (at the poles gradients will be 'nan').

From this, the following ancillaries generated:

 - filtered mean (m01s00i033)
 - standard deviation (m01s00i034)
 - sigma_xx (m01s00i035)
 - sigma_yy (m01s00i036)
 - sigma_xy (m01s00i037)

Additionally the following LM97 parameters are generated using the fields
above:

 - Orographic sub-grid anisotropy (m01s06i248)
 - Orographic sub-grid orientation (m01s06i249)
 - Orographic sub-grid slope (m01s06i250)

See the following source on how the LM97 parameters are derived:

https://code.metoffice.gov.uk/doc/um/vn11.0/papers/umdp_022.pdf#subsection.3.1

"""
import functools
import logging

import ants
import ants.decomposition as decomp
import ants.io.save as save
import iris
import iris.analysis.calculus
import numpy as np
import orography_utils

iris.FUTURE.save_split_attrs = True
_LOGGER = logging.getLogger(__name__)


def calc_grad_xy(grad_x, grad_y):
    grad_xy = iris.analysis.maths.multiply(grad_x, grad_y)
    grad_xy.long_name = "Orographic Gradient XY Component"
    grad_xy.units = "m2/m2"
    grad_xy.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i036")
    return grad_xy


def calc_grad_yy(grad_y):
    grad_yy = grad_y ** 2
    grad_yy.long_name = "Orographic Gradient YY Component"
    grad_yy.units = "m2/m2"
    grad_yy.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i037")
    return grad_yy


def calc_grad_xx(grad_x):
    grad_xx = grad_x ** 2
    grad_xx.long_name = "Orographic Gradient XX Component"
    grad_xx.units = "m2/m2"
    grad_xx.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i035")
    return grad_xx


def stdev(source, src_mean, regridder):
    """
    Calculate the standard deviation of the source within the target grid box.

    """
    source **= 2
    src_mean **= 2
    awm = regridder(source)
    awm -= src_mean
    # Account for precision issues resulting in negative values.
    awm.data = np.abs(awm.data)
    awm **= 1 / 2.0
    return awm


def calc_lm97_params(sigma_xx, sigma_yy, sigma_xy):
    """
    Derive the LM97 sub-grid scale parameters.

    Parameters
    ----------
    sigma_xx : :class:`~iris.cube.Cube`
        Orographic Gradient XX Component.
    sigma_yy : :class:`~iris.cube.Cube`
        Orographic Gradient YY Component
    sigma_xy : :class:`~iris.cube.Cube`
        Orographic Gradient XY Component

    Returns
    -------
    : :class:`~iris.cube.Cube`, :class:`~iris.cube.Cube`,
            :class:`~iris.cube.Cube`
        - Orographic sub-grid anisotropy (m01s06i248)
        - Orographic sub-grid orientation (m01s06i249)
        - Orographic sub-grid slope (m01s06i250)

    See Also
    --------
    See the following technical document which describes the above returns:
    https://code.metoffice.gov.uk/doc/um/vn11.0/papers/umdp_022.pdf#\
    subsection.3.1

    """
    smallp = np.finfo("float32").eps / 100.0

    tmp = sigma_xx.data + sigma_yy.data
    tmp1 = sigma_xx.data - sigma_yy.data
    tmp2 = (tmp1 ** 2 + (4 * sigma_xy.data ** 2)) ** (1 / 2.0)
    sso_anisotropy_numer = tmp - tmp2
    sso_anisotropy_denom = tmp + tmp2
    sso_anisotropy = sso_anisotropy_numer / sso_anisotropy_denom
    sso_anisotropy **= 1 / 2.0
    sso_anisotropy = sigma_xx.copy(sso_anisotropy)
    sso_anisotropy.rename("Orographic sub-grid anisotropy")
    sso_anisotropy.data[sso_anisotropy_denom <= smallp] = 1
    sso_anisotropy.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi(
        "m01s06i248"
    )

    psi = sigma_xy.copy()
    psi.data = 0.5 * np.arctan2(2 * sigma_xy.data, tmp1)
    psi.rename("Orographic sub-grid orientation")
    psi.data[(np.abs(sigma_xy.data) <= smallp) * (np.abs(tmp1) <= smallp)] = 0
    psi.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s06i249")

    mean_slope = (0.5 * sso_anisotropy_denom) ** (1 / 2.0)
    mean_slope = sigma_xx.copy(mean_slope)
    mean_slope.rename("Orographic sub-grid slope")
    mean_slope.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s06i250")

    return sso_anisotropy, psi, mean_slope


def calc_fields(mean, polar_cutoff, source, ocean_mask):
    """
     Calculate the sub grid scale orography (SSO) fields.

    Parameters
    ----------
    polar_cutoff : float
        Absolute latitude value in degrees beyond which to set grad_x and
        grad_y to 0.
    source : :class:`iris.cube.Cube`
        Orography source dataset.
    ocean_mask : :class:`iris.cube.Cube`
        Target landsea mask.

    Returns
    -------
    : :class:`iris.cube.Cube`, :class:`iris.cube.Cube`, :class:`iris.cube.Cube`
            :class:`iris.cube.Cube`, :class:`iris.cube.Cube`,
            :class:`iris.cube.Cube`, :class:`iris.cube.Cube`,
            :class:`iris.cube.Cube`
        filtered mean (m01s00i033), standard deviation (m01s00i034),
        sigma_xx (m01s00i035), sigma_yy (m01s00i036), sigma_xy (m01s00i037),
        anisotropy of the SSO, angle between the low-level wind and the major
        axis of the topography, the mean slope of the SSO along the major axis.

    See Also
    --------
    https://code.metoffice.gov.uk/doc/um/vn11.0/papers/umdp_015.pdf#\
    subsubsection.2.1.3

    """

    tymin, tymax = (
        ocean_mask.coord(axis="y").points.min(),
        ocean_mask.coord(axis="y").points.max(),
    )
    txmin, txmax = (
        ocean_mask.coord(axis="x").points.min(),
        ocean_mask.coord(axis="x").points.max(),
    )
    slices = ants.utils.cube.get_slices(
        mean, [tymin, tymax], [txmin, txmax], pad_width=0
    )
    assert len(slices) == 1, "mean isn't on the same grid as the ocean_mask"
    mean = mean[slices[0]]
    assert (mean.coord(axis="x").points == ocean_mask.coord(axis="x").points).all()
    print(f"mean y points: {mean.coord(axis='y').points}")
    print(f"mask y points: {ocean_mask.coord(axis='y').points}")
    assert (mean.coord(axis="y").points == ocean_mask.coord(axis="y").points).all()

    grad_x, grad_y = ants.analysis.calc_grad(source)
    lats = source.coord(axis="y").points

    if polar_cutoff is not None:
        pole_cut_off = polar_cutoff + 1e-10
        grad_x.data[np.abs(lats) > pole_cut_off, :] = 0.0
        grad_y.data[np.abs(lats) > pole_cut_off, :] = 0.0

    grad_xx = calc_grad_xx(grad_x)
    grad_yy = calc_grad_yy(grad_y)
    grad_xy = calc_grad_xy(grad_x, grad_y)

    regridder = ants.regrid.GeneralRegridScheme(horizontal_scheme="TwoStage").regridder(
        source, ocean_mask
    )
    sigma_xx = regridder(grad_xx)
    sigma_yy = regridder(grad_yy)
    sigma_xy = regridder(grad_xy)

    # Calculate the (orog_height - mean gradient field) weighted average.
    source_smg_tar = regridder(source)
    stddev = stdev(source, source_smg_tar, regridder)
    stddev.rename("orography standard deviation")
    stddev.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i034")

    # Ensure consistency with the supplied ocean mask.
    for cube in [mean, stddev, sigma_xx, sigma_yy, sigma_xy]:
        orography_utils.make_consistent_with_mask(cube, ocean_mask)

    mean.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i033")

    # Set nan values to 0 (this will be at the poles...).  We could constrain
    # this to the pole more explicitly perhaps.
    for cube in [stddev, sigma_xx, sigma_yy, sigma_xy]:
        cube.data[np.isnan(cube.data)] = 0

    sso_anisotropy, psi, mean_slope = calc_lm97_params(sigma_xx, sigma_yy, sigma_xy)
    return (mean, stddev, sigma_xx, sigma_yy, sigma_xy, sso_anisotropy, psi, mean_slope)


def calc_source_sub_mean_grad(mean, source):
    # Calculate an estimate of the mean gradient field by performing a bilinear
    # interpolation of the mean back onto the source grid.
    # This is its own separate decomposed operation as we want to iterate over
    # the source and fetch the corresponding mean chunk to ensure that we
    # generate a target with grid definition identical to the source.
    regridder = ants.regrid.rectilinear.Linear(extrapolation_mode="linear")
    mean_src = mean.regrid(source, regridder)
    mean_src.units = "m"
    source.units = "m"
    return source - mean_src


def load_data(filt_source, target, mean):
    src_cube_filt = ants.io.load.load_cube(filt_source, "surface_altitude_filtered")
    tgt_cube = ants.io.load.load_cube(target, "land_binary_mask")
    mean_cube = ants.io.load.load_cube(mean, iris.AttributeConstraint(STASH="m01s00i033"))
    mean_cube.rename("orography mean height")
    ants.utils.cube.guess_horizontal_bounds(src_cube_filt)
    ants.utils.cube.guess_horizontal_bounds(tgt_cube)
    ants.utils.cube.guess_horizontal_bounds(mean_cube)
    # Remove any grid staggering.
    src_cube_filt.attributes.clear()
    tgt_cube.attributes.clear()
    mean_cube.attributes.clear()

    mean_cube.units = "m"
    mean_cube.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i033")

    # Subset for regional grids
    if not ants.utils.cube.is_global(tgt_cube):
        src_cube_filt = src_cube_filt.extract(ants.ExtractConstraint(tgt_cube))

    return src_cube_filt, tgt_cube, mean_cube


def main(filt_source, target, meanname, output, mean_orography_output, netcdf_only, polar_cutoff):

    orog_height_filt, target, mean = load_data(filt_source, target, meanname)

    ocean = orography_utils.gen_ocean_mask(target)

    # estimate (orog_height_filt - mean gradient field)
    source_smg = decomp.decompose(calc_source_sub_mean_grad, mean, orog_height_filt)

    operation = functools.partial(calc_fields, mean, polar_cutoff)

    results = decomp.decompose(operation, source_smg, ocean)

    save.netcdf(results, output)
    if not netcdf_only:
        save.ancil(results, output)

    # separately save mean after ocean mask
    orography_utils.make_consistent_with_mask(mean, ocean)

    save.netcdf(mean, mean_orography_output)
    if not netcdf_only:
        save.ancil(mean, mean_orography_output)


if __name__ == "__main__":
    parser = ants.AntsArgParser(target_lsm=True)
    parser.add_argument(
        "--mean-orography-input",
        required=True,
        type=ants.config.filepath_readable,
        help=("Mean orography input file."),
    )

    parser.add_argument(
        "--mean-orography-output",
        type=ants.config.dirpath_writeable,
        required=True,
        help=("File to save mean orography ancil to."),
    )

    parser.add_argument(
        "--polar-cutoff",
        nargs="?",
        const=88.0,
        type=float,
        help=(
            "Apply a polar cutoff above to grad_x and grad_y calculations "
            "by zeroing values above whose absolute latitude value is above "
            "a cutoff. If a cutoff value is not provided then a default "
            "cutoff threshold of 88 degrees will be used."
        ),
    )

    args = parser.parse_args()

    if args.output == args.mean_orography_output:
        raise parser.error(
            "Mean orography must be saved to a different path than "
            "the other ancillaries."
        )

    main(
        args.sources,
        args.target_lsm,
        args.mean_orography_input,
        args.output,
        args.mean_orography_output,
        args.netcdf_only,
        args.polar_cutoff,
    )
