# ancil_clim_dms

The `ancil_clim_dms` app produces the dimethyl sulphide (DMS) sea-water concentration ancillary used by the CLASSIC aerosol scheme. It takes the Lana et al. ocean DMS climatology and regrids it onto the target land-sea mask grid, then makes it consistent with that mask so that only ocean points retain data. Any remaining masked (land) points are replaced with 0, since there is no DMS concentration over land.

## Env Arguments

#### Inputs
* `DMS_SOURCE`: Source DMS climatology. Must be a NetCDF dataset containing a `mole_concentration_of_dimethyl_sulfide_in_sea_water` variable.
* `TARGET_LSM`: Target land-sea mask to regrid onto and make the DMS field consistent with. Must be a NetCDF dataset containing a `land_binary_mask` variable.

#### Outputs
* `OUTPUT`: Path to write the output to, in UM PP format and also in NetCDF format.

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
