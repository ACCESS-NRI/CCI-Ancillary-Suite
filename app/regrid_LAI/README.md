# regrid_LAI

The `regrid_LAI` app regrids the preprocessed MODIS-derived total (all-vegetation) leaf area index (LAI) climatology onto the target land-sea mask grid, in preparation for use by `ancil_LAI`, which derives the per-plant-functional-type LAI from this total LAI field. It invokes `ancil_general_regrid.py`, a generic ANTS regridding wrapper shared by several tasks in this suite, which regrids source data onto a target grid using `ants.regrid.GeneralRegridScheme`.

Because a target land-sea mask is supplied, the regridded result is made consistent with that mask: any masked/missing points are filled with valid neighbouring values via a nearest-neighbour (spiral) search, so that the output honours the coastlines of the target grid. The `--netcdf-only` option is used, so the output is only written in NetCDF format.

## Env Arguments

#### Inputs
* `LAI_SOURCE`: Preprocessed MODIS total LAI source dataset.
* `TARGET_LSM`: Target land-sea mask defining the destination grid, from `ancil_lct`.

#### Outputs
* `OUTPUT`: Path to write the regridded LAI to, in NetCDF format only. Consumed as `LAI_SOURCE` by `ancil_LAI`.

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
