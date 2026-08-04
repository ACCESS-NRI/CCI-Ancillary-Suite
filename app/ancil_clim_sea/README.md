# ancil_clim_sea

The `ancil_clim_sea` app regrids the GlobColour ocean colour ("sea colour") climatology onto the target land-sea mask grid to produce the `qrclim.sea` ancillary. It invokes `ancil_general_regrid.py`, a generic ANTS regridding wrapper shared by several tasks in this suite, which regrids source data onto a target grid using `ants.regrid.GeneralRegridScheme`.

Because a target land-sea mask is supplied, the regridded result is made consistent with that mask: any masked/missing points are filled with valid neighbouring values via a nearest-neighbour (spiral) search, so that the output honours the coastlines of the target grid.

## Env Arguments

#### Inputs
* `SEA_COLOUR_SOURCE`: Source ocean colour climatology, from the GlobColour dataset.
* `TARGET_LSM`: Target land-sea mask defining the destination grid, from `ancil_lct`.

#### Outputs
* `OUTPUT`: Path to write the output to, in UM PP format and also in NetCDF format.

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
