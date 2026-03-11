#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
"""
Dimethylsulfide (DMS) concentration generation
**********************************************
- Perform regrid to the specified landsea mask grid.
- Make consistent with the provided mask.
- Replace mask with 0, ensuring values of 0 DNS over land water (where there
  is no DNS).

"""
import ants
import ants.decomposition as decomp
import ants.io.save as save
import numpy as np


def load(sources, target_lsm):
    lsm_cube = ants.io.load.load_cube(target_lsm, "land_binary_mask")
    dms_cube = ants.io.load.load_cube(
        sources, "mole_concentration_of_dimethyl_sulfide_in_sea_water"
    )
    return lsm_cube, dms_cube


def operation(source, target):
    scheme = ants.regrid.GeneralRegridScheme(horizontal_scheme="TwoStage")
    result = source.regrid(target, scheme)
    return result


def main(sources, output, target_lsm, netcdf_only):
    lsm_cube, dms_cube = load(sources, target_lsm)
    res_cube = decomp.decompose(operation, dms_cube, lsm_cube)
    ants.analysis.make_consistent_with_lsm(res_cube, lsm_cube, False)
    if np.ma.isMaskedArray(res_cube.data):
        res_cube.data = res_cube.data.filled(0)

    save.netcdf(res_cube, output)
    if not netcdf_only:
        save.ancil(res_cube, output)

    return res_cube


if __name__ == "__main__":
    parser = ants.AntsArgParser(target_lsm=True)
    args = parser.parse_args()
    main(
        args.sources,
        args.output,
        args.target_lsm,
        args.netcdf_only,
    )
