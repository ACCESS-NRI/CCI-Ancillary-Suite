# regrid_dust_aerosol

The `regrid_dust_aerosol` app regrids the GA6.0 mineral dust aerosol climatology onto the target model grid, producing the `qrclim.dust` ancillary. It is one of seven per-species aerosol regridding apps (`regrid_{biog,biom,blck,dust,ocff,seasalt,sulp}_aerosol`) that run when `AEROSOLS="PP"` in `rose-suite.conf`; when `AEROSOLS="NetCDF"`, the single `regrid_aerosol` app is run instead. It invokes `ancil_general_regrid.py`, a generic ANTS regridding wrapper shared by several tasks in this suite, which regrids source data onto a target grid using `ants.regrid.GeneralRegridScheme`.

Unlike the tasks that regrid onto a land-sea mask, this app supplies a target grid made up of two files: a horizontal grid definition (`TARGET_OROGRAPHIC_LSM`, produced by `ancil_orographic_wavedrag`) and a vertical levels namelist (`VERTICAL_DISCRETIZATION`), so both horizontal and vertical regridding are performed, each using a linear scheme. The source climatology is read directly from the ancillary master directory, so no pre-processing task is required. The result is saved as a UM ancillary file, alongside a NetCDF copy.

## Env Arguments

#### Inputs
* `DUST_SOURCE`: Source mineral dust aerosol climatology, from `${ANCIL_MASTER}/atmos/master/aero_clims/GA6.0_antie/v1/qrclim.dust.nc`.
* `TARGET_OROGRAPHIC_LSM`: Horizontal component of the target grid, from `ancil_orographic_wavedrag`. Set in `flow.cylc`.
* `VERTICAL_DISCRETIZATION`: Vertical component of the target grid. Should be a namelist in the UM vertical discretization format. Set globally in `rose-suite.conf`.

#### Outputs
* `OUTPUT`: Path to write the regridded mineral dust aerosol climatology to, as a UM ancillary (via `ants.io.save.ancil`) plus a NetCDF copy (via `ants.io.save.netcdf`).

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
