# regrid_seaice_reynolds

The `regrid_seaice_reynolds` app regrids the Reynolds sea ice concentration archive (spanning September 1981 to December 2019) onto the target ocean grid, producing the `qrclim.seaice` ancillary. It invokes `ancil_general_regrid.py`, a generic ANTS regridding wrapper shared by several tasks in this suite, which regrids source data onto a target grid using `ants.regrid.GeneralRegridScheme`.

Source data outside the `BEGIN_YEAR`-`END_YEAR` range is discarded before regridding. Because a target ocean mask is supplied, the regridded result is made consistent with that mask: any masked/missing points are filled with valid neighbouring values via a nearest-neighbour (spiral) search, so that the output honours the coastlines of the target grid.

## Env Arguments

#### Inputs
* `SEAICE_SOURCE`: Source Reynolds sea ice concentration dataset. The default path includes the `CALENDAR` setting (from `rose-suite.conf`) to select the calendar-specific archive.
* `TARGET_OCEAN`: Target ocean mask defining the destination grid, from `ancil_lct`.
* `BEGIN_YEAR`: Start year (inclusive) of the source data range to include; data prior to this year is discarded. Set globally in `rose-suite.conf`.
* `END_YEAR`: End year (inclusive) of the source data range to include; data after this year is discarded. Set globally in `rose-suite.conf`.

#### Outputs
* `OUTPUT`: Path to write the output to, in UM PP format and also in NetCDF format.

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
