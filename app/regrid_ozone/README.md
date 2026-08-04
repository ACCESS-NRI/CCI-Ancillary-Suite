# regrid_ozone

The `regrid_ozone` app regrids the CMIP5 SPARCex ozone climatology (1994-2005) onto the target model grid, producing the `qrclim.ozone` ancillary. It invokes `ancil_general_regrid.py`, a generic ANTS regridding wrapper shared by several tasks in this suite, which regrids source data onto a target grid using `ants.regrid.GeneralRegridScheme`.

As with `regrid_aerosol`, the target grid is made up of two files: a horizontal grid definition (`TARGET_OROGRAPHIC_LSM`, produced by `ancil_orographic_wavedrag`) and a vertical levels namelist (`VERTICAL_DISCRETIZATION`), so both horizontal and vertical regridding are performed.

## Env Arguments

#### Inputs
* `OZONE_SOURCE`: Source ozone climatology, from the CMIP5 SPARCex dataset.
* `TARGET_OROGRAPHIC_LSM`: Horizontal component of the target grid, from `ancil_orographic_wavedrag`.
* `VERTICAL_DISCRETIZATION`: Vertical component of the target grid. Should be a namelist in the UM vertical discretization format. Set globally in `rose-suite.conf`.

#### Outputs
* `OUTPUT`: Path to write the output to, in UM PP format and also in NetCDF format.

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
