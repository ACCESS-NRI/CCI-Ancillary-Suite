#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
Land cover type fraction post-processor for C4 Still dataset injection
**********************************************************************
- Use the provided C4 Still dataset into the land cover type fraction datase
  by the following relation::

        vegfrac_C4_grass = min(C4_still, vegfrac_C3_grass)
        vegfrac_C3_grass = vegfrac_C3_grass - vegfrac_C4_grass

    where the input `vegfrac_C3_grass` represents the combined C3 and C4
    grass fraction to be split.

"""
import ants
import ants.decomposition as decomp
import ants.io.save as save
import iris
import numpy as np


def load(lct, c4_source):
    c4_cube = ants.io.load.load_cube(c4_source)
    stash_con = iris.AttributeConstraint(STASH="m01s00i216")
    lct_cube = ants.io.load.load_cube(lct, constraint=stash_con)
    return lct_cube, c4_cube


def derive_c4_contributing(c4_source, lct):
    """
    Inject C4 data into the land cover type fraction dataset.

    The C4 is incorporated as follows::

        lct_C4_grass = min(C4_still, lct_C3_grass)
        lct_C3_grass = lct_C3_grass - lct_C4_grass

    where the input `lct_C3_grass` represents the combined C3 and C4
    grass fraction to be split.

    Parameters
    ----------
    lct : :class:`~iris.cube.Cube`
        Vegetation fraction dataset.
    c4_source : :class:`~iris.cube.Cube`
        c4 Still etc al. source data used to distinguish between C3/C4 JULES
        classes.
    c3level : int
        Index of C3 level in the land cover type fraction dataset (default 3)
    c4level : int
        Index of C4 level in the land cover type fraction dataset (default 4)

    Returns
    -------
    : :class:`~iris.cube.Cube`
        Vegetation fraction with C3/C4 distinction.

    """
    c4_source = ants.analysis.mean(c4_source, lct)
    # Convert percentage to fraction
    c4_source.data = c4_source.data / 100.0

    # Make sure the C4 still has data consistent with the CCI mask (ocean).
    mask = lct[0].copy(np.ma.getmaskarray(lct[0].data))
    filler = ants.analysis.FillMissingPoints(c4_source, target_mask=mask)
    filler(c4_source)

    c3_slice = ants.analysis.cover_mapping.fetch_lct_slices(lct, args.c3level)
    c4_slice = ants.analysis.cover_mapping.fetch_lct_slices(lct, args.c4level)

    # Check there isn't any pre-existing C4 grass
    if lct.data[c4_slice].any():
        msg = (
            "There appears to be some pre-existing C4 grass fraction "
            "present, perhaps you don't need to inject the C4 grass from "
            "the Still or perhaps you have done so already?"
        )
        raise ValueError(msg)

    # Handle mask difference with lct
    lct.data[c4_slice] = np.minimum(c4_source.data, lct.data[c3_slice])
    lct.data[c3_slice] = lct.data[c3_slice] - lct.data[c4_slice]
    lct.rename("vegetation_area_fraction")
    return lct


def main(lct, c4_source, output, netcdf_only):
    lct, c4_source = load(lct, c4_source)
    lct = decomp.decompose(derive_c4_contributing, c4_source, lct)

    if not netcdf_only:
        save.ancil(lct, output)
    save.netcdf(lct, output)
    return lct


def _get_parser():
    parser = ants.AntsArgParser()
    c4_help = (
        "Pathname of ISLSCP II C4 Vegetation Percentage source to "
        "derive C3/C4 distinction. See "
        "https://daac.ornl.gov/ISLSCP_II/guides/c4_percent_1deg.html"
    )
    parser.add_argument("--c4-source", type=str, help=c4_help, required=True)
    parser.add_argument("--c3level", type=int, help="index of C3 level", default=3)
    parser.add_argument("--c4level", type=int, help="index of C4 level", default=4)
    return parser


if __name__ == "__main__":
    args = _get_parser().parse_args()
    main(args.sources, args.c4_source, args.output, args.netcdf_only)
