#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
"""
- The mean orography is derived:
    - Conservatively regridding the unfiltered source orographic height field
      onto the target grid.
    If not using blended orography:
    - Filter the mean using the Raymond filter with epsilon defined by the
      command line argument --epsilon. The epsilon value is defined at the
      equator (isotropy is maintained where possible, see
      :func:`ants.analysis.filters.raymond`). For global use, the default
      value (0.01) should be used. For regional use, this parameter may need
      to be tuned for model stability.
    If using blended orography:
    If computing using max_slope_from_smooth:
    - Get min and max epsilon values and blend them using blending script
    If computing using max_slope:
    - Loop through defined epsilon values until you find the one that has a maximum
      slope < the given max_slope
      If this condition is not satisfied by any, use the maximum epsilon value and
      only perform looping over max_smooth iterations
    - Blend the orography over regions where the slope is larger than the
      max_slope value
"""

import logging

import ants
import ants.decomposition as decomp
import ants.io.save as save
import iris
import iris.analysis.calculus
import numpy as np
import orography_utils
from ants.analysis.filters import raymond

_LOGGER = logging.getLogger(__name__)


def filter_mean(mean, epsilon_value=0.01):
    """
    Filter mean using the raymond filter
    We cannot use this filter on decomposed pieces,
    so this function must be called before decomposition
    """
    tgt_crs = mean.coord(axis="x").coord_system.as_ants_crs()
    isotropic = (tgt_crs == ants.coord_systems.UM_SPHERE.crs) and (
        ants.utils.cube.is_global(mean)
    )
    raymond(mean, epsilon=epsilon_value, isotropic=isotropic)
    mean.data[mean.data < 0] = 0
    return ants.utils.cube.defer_cube(mean)


def load_data(source, target):
    src_cube = ants.io.load.load_cube(source)
    tgt_cube = ants.io.load.load_cube(target, "land_binary_mask")
    ants.utils.cube.guess_horizontal_bounds(src_cube)
    ants.utils.cube.guess_horizontal_bounds(tgt_cube)
    # Remove any grid staggering.
    src_cube.attributes.clear()
    tgt_cube.attributes.clear()
    src_cube.units = "m"
    return src_cube, tgt_cube


def find_smoothest_orography(mean, max_slope, epsilon, max_slope_from_smooth):
    """
    Generate uniformly smoothed orography dataset compliant with specified
    maximum slope, for use in blend for targeted smoothing.
    """
    for eps in sorted(epsilon):
        print("Trying with epsilon =", eps)
        smooth_orography = filter_mean(mean, eps)
        smooth_grad = orography_utils.get_total_gradient(smooth_orography)
        print(
            "smooth orography gradient, epsilon = %s, gradient = %s"
            % (eps, smooth_grad.max())
        )

        isteep = np.where(smooth_grad > max_slope)
        if isteep[0].size == 0:
            break
        else:
            print("smoothest orography exceeds max_slope, using maximum epsilon value")
            # and setting max_slope_from_smooth = True
            # max_slope_from_smooth = True
    return smooth_orography, eps


def main(
    sources,
    target,
    output,
    epsilon,
    blend_orography,
    max_slope,
    max_slope_from_smooth,
    max_smooth,
    dilation,
):

    orog_height, target = load_data(sources, target)

    mean = decomp.decompose(orography_utils.calc_mean, orog_height, target)

    if blend_orography:
        fine_orography = filter_mean(mean, min(epsilon))
        if max_slope_from_smooth:
            # Compute the blended orography using the maximum slope from
            # the specified maximum epsilon value
            smooth_orography = filter_mean(mean, max(epsilon))
            mean = orography_utils.blended_orography(
                fine_orography,
                smooth_orography,
                max_slope,
                max_slope_from_smooth,
                max_smooth,
                dilation=dilation,
            )
        else:
            # Compute using max_slope
            # First check that max slope not already satisfied
            fine_grad = orography_utils.get_total_gradient(fine_orography)
            print("fine orography gradient = ", fine_grad.max())
            isteep = np.where(fine_grad > max_slope)
            if isteep[0].size != 0:
                # Find a suitable epsilon value for smoothest orography
                smooth_orography, smooth_epsilon = find_smoothest_orography(
                    mean, max_slope, epsilon, max_slope_from_smooth
                )
                mean = orography_utils.blended_orography(
                    fine_orography,
                    smooth_orography,
                    max_slope,
                    max_slope_from_smooth,
                    max_smooth,
                    dilation=dilation,
                )
            else:
                mean = fine_orography
                print("max_slope not exceeded, using minimum epsilon value")
    else:
        mean = filter_mean(mean, min(epsilon))

    # save pre-ocean masking for file to be supplied to wavedrag code
    mean.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i033")
    save.netcdf(mean, output)
    save.ancil(mean, output)


if __name__ == "__main__":
    parser = ants.AntsArgParser(target_lsm=True)
    parser.add_argument(
        "--blend-orography",
        action="store_true",
        help=(
            "Logical to generate the mean orography by blending orographies."
            "Defaults to False."
        ),
    )
    parser.add_argument(
        "--dilation",
        action="store_true",
        help=(
            "Logical to use dilation."
            "Defaults to False."
        ),
    )
    default_epsilon = [0.01, 0.05, 0.1, 1.0, 2.0, 5.0, 8.0, 10.0]
    parser.add_argument(
        "--epsilon",
        nargs="*",
        type=float,
        help=(
            "Epsilon value(s) for the Raymond filter, defined at the equator."
            f"Defaults to {default_epsilon}."
        ),
        default=default_epsilon,
    )

    # Default max_slope value is chosen to be the maximum slope of the N640
    # mean orography that has been filtered with epsilon=10.0
    default_max_slope = 0.02
    parser.add_argument(
        "--max-slope",
        type=float,
        help=("Value of the maximum slope" f"Defaults to {default_max_slope}."),
        default=default_max_slope,
    )
    # Default max_smooth value was found by trial to be ample for convergence in
    # the cases where the algorithm has been applied to date.
    default_max_smooth = 30
    parser.add_argument(
        "--max-smooth",
        type=float,
        help=(
            "Value of the maximum number smoothing iterations"
            f"Defaults to {default_max_smooth}."
        ),
        default=default_max_smooth,
    )
    parser.add_argument(
        "--max-slope-from-smooth",
        action="store_true",
        help=("Logical to use the maximum slope from the smoothed orography"),
    )
    args = parser.parse_args()
    if all([args.blend_orography, len(args.epsilon) < 2]):
        raise RuntimeError(
            "When --blend_orography option is selected, at least two epsilon"
            "values must be entered in order of smallest (finest) to largest"
            "(coarsest)."
        )
    if all(
        [args.blend_orography, not args.max_slope_from_smooth, len(args.epsilon) == 2]
    ):
        raise RuntimeWarning(
            "When --blend_orography and --max_slope is specified,"
            "more than 2 epsilon values should be listed, otherwise method is"
            "identical to selecting --max-slope-from-smooth. "
            "Therefore, setting max_slope_from_smooth = True."
        )
        max_slope_from_smooth = True
    main(
        args.sources,
        args.target_lsm,
        args.output,
        args.epsilon,
        args.blend_orography,
        args.max_slope,
        args.max_slope_from_smooth,
        args.max_smooth,
        args.dilation,
    )
