# ancil_river_storage

The `ancil_river_storage` app derives the river storage ancillary (`qrclim.rivstor`) by combining the river storage source data with the river direction field produced by `ancil_river_routing`. The precise steps performed by the task are:

1. The river storage source is regridded onto the same grid as the river direction ancillary, using an area-weighted mean.
2. The mask of the direction cube (i.e. points with no defined river direction) is applied to the regridded storage cube, so that storage values are only retained where a river direction is defined.

## Env Arguments

#### Inputs
* `ROUTING_SOURCE`: River routing ancillary (`qrparm.rivseq`) produced by `ancil_river_routing`. The river direction field (STASH `m01s00i152`) is extracted from this source and used both to define the target regrid and to mask the output.
* `STORAGE_SOURCE`: River storage source data produced by `ancil_river_storage_preproc`, containing the storage field with STASH `m01s00i153`.

#### Outputs
* `OUTPUT`: Path to write the masked, regridded river storage field to, in UM PP format, plus NetCDF format with the `.nc` suffix added.

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
