# ancil_river_storage_preproc

The `ancil_river_storage_preproc` app corrects metadata on the source Fekete river storage dataset so that it is interpretable as a periodic monthly-mean time series in the final ancillary. Specifically, the source PP fields (STASH `m00s26i001`) are re-tagged as monthly means representative of the 1950-2000 period, with the time metadata (`lbtim`, `lbdatd`, `lbyr`/`lbyrd`, `lbmond`/`lbmind`/`lbsec` etc.) rewritten so that a full 12-month source is recognised as periodic. If the source data has already been corrected (or never needed correcting), the correction is skipped and a warning is raised instead. The resulting cube's STASH is re-tagged to `m01s00i153` (river routing storage) and its time coordinate is annotated with a `representative_period` attribute of "1950 - 2000".

This app feeds its NetCDF output into `ancil_river_storage` as `STORAGE_SOURCE`. A comment in the app's `rose-app.conf` notes that this step appears independent of resolution or model configuration, and questions whether it needs to be run as part of the suite at all rather than once upstream.

## Env Arguments

#### Inputs
* `STORAGE_SOURCE`: Source river storage data (Fekete dataset) in UM PP format, containing the river storage field with STASH `m00s26i001`.

#### Outputs
* `OUTPUT`: Path to write the corrected river storage field to. Only written in NetCDF format (the `--netcdf-only` flag is set in the command), for consumption by `ancil_river_storage`.

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
