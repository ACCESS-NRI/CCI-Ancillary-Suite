# ancil_river_routing_preproc

The `ancil_river_routing_preproc` app corrects metadata on the source TRIP river routing dataset (river routing sequence and direction fields) so that it can be safely interpreted downstream. Specifically, it modifies the source time points so that they do not exceed 30 days in a month, which would otherwise be invalid for the model's 360-day calendar. If the source data has already been corrected (or never needed correcting), the correction is skipped and a warning is raised instead.

This app feeds its NetCDF output into `ancil_river_routing` as `ROUTING_SOURCE`. A comment in the app's `rose-app.conf` notes that this step appears independent of resolution or model configuration, and questions whether it needs to be run as part of the suite at all rather than once upstream.

## Env Arguments

#### Inputs
* `ROUTING_SOURCE`: Source river routing data in UM PP format, containing the river routing sequence (STASH `m01s00i151`) and river routing direction (STASH `m01s00i152`) fields.

#### Outputs
* `OUTPUT`: Path to write the corrected river routing sequence and direction fields to. Only written in NetCDF format (the `--netcdf-only` flag is set in the command), for consumption by `ancil_river_routing`.

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
