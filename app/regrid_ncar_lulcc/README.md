# regrid_ncar_lulcc

The `regrid_ncar_lulcc` app regrids the NCAR Land Use and Land Cover Change (LULCC) plant functional type (PFT) fraction dataset onto the target land-sea mask grid, in preparation for use by `ancil_split_grass_types`, which uses it to split the short vegetation classes of the CCI land cover into more detailed sub-types (C3, C4, arctic grasses and shrub types).

Before regridding, `ncks` and `ncrename` are used to extract the `LON`, `LAT` and `PCT_PFT` variables from the source file and rename them to `lon`, `lat` and `vegetation_area_fraction` respectively, so that the data is recognised by ANTS. The extracted data is then regridded using `ancil_general_regrid.py`, a generic ANTS regridding wrapper shared by several tasks in this suite, which regrids source data onto a target grid using `ants.regrid.GeneralRegridScheme`. This app uses a variant of the script that additionally assumes a default (`GeogCS`) coordinate system for any source or target coordinates that are missing one, and performs the regrid directly rather than via domain decomposition.

Because a target land-sea mask is supplied, the regridded result is made consistent with that mask: any masked/missing points are filled with valid neighbouring values via a nearest-neighbour (spiral) search, so that the output honours the coastlines of the target grid. The `--netcdf-only` option is used, so the output is only written in NetCDF format.

## Env Arguments

#### Inputs
* `LULCC_SOURCE`: Source NCAR LULCC dataset, containing `LON`, `LAT` and `PCT_PFT` variables.
* `TARGET_LSM`: Target land-sea mask defining the destination grid, from `ancil_lct`.

#### Outputs
* `OUTPUT`: Path to write the regridded PFT fractions to, in NetCDF format only. Consumed as `GRASS_SOURCE` by `ancil_split_grass_types`.

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
