# regrid_C4_fraction

The `regrid_C4_fraction` app regrids the CCI C4 grass percentage dataset (derived from the ISLSCP II C4 Vegetation Percentage, Still et al.) onto the target land-sea mask grid, in preparation for use by `ancil_C4_fraction`, which uses it to split the combined C3/C4 grass fraction in the land cover type ancillary. It invokes `ancil_general_regrid.py`, a generic ANTS regridding wrapper shared by several tasks in this suite, which regrids source data onto a target grid using `ants.regrid.GeneralRegridScheme`.

Because a target land-sea mask is supplied, the regridded result is made consistent with that mask: any masked/missing points are filled with valid neighbouring values via a nearest-neighbour (spiral) search, so that the output honours the coastlines of the target grid. The `--netcdf-only` option is used, so the output is only written in NetCDF format.

## Env Arguments

#### Inputs
* `C4_SOURCE`: Source CCI C4 grass percentage dataset.
* `TARGET_LSM`: Target land-sea mask defining the destination grid, from `ancil_lct`.

#### Outputs
* `OUTPUT`: Path to write the regridded C4 fraction to, in NetCDF format only. Consumed as `C4_SOURCE` by `ancil_C4_fraction`.

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
