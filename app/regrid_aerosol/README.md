# regrid_aerosol

The `regrid_aerosol` app regrids the pre-processed GLOMAP MODE aerosol climatology (produced by `ancil_aeroclim_preproc`) onto the target model grid, producing the `qrclim.aerosols` ancillary used by the UKCA aerosol scheme. It invokes `ancil_general_regrid.py`, a generic ANTS regridding wrapper shared by several tasks in this suite, which regrids source data onto a target grid using `ants.regrid.GeneralRegridScheme`.

Unlike the tasks that regrid onto a land-sea mask, this app supplies a target grid made up of two files: a horizontal grid definition (`TARGET_OROGRAPHIC_LSM`, produced by `ancil_orographic_wavedrag`) and a vertical levels namelist (`VERTICAL_DISCRETIZATION`), so both horizontal and vertical regridding are performed. The result is saved using the UKCA-specific NetCDF saver (`--save-ukca`), and a subsequent `ncatted` call sets the `update_type` global attribute on the output to `2`.

## Env Arguments

#### Inputs
* `AEROSOL_SOURCE`: Pre-processed aerosol climatology, from `ancil_aeroclim_preproc`.
* `TARGET_OROGRAPHIC_LSM`: Horizontal component of the target grid, from `ancil_orographic_wavedrag`.
* `VERTICAL_DISCRETIZATION`: Vertical component of the target grid. Should be a namelist in the UM vertical discretization format. Set globally in `rose-suite.conf`.

#### Outputs
* `OUTPUT`: Path to write the regridded aerosol climatology to, in UKCA-specific NetCDF format (via `ants.io.save.ukca_netcdf`).

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
