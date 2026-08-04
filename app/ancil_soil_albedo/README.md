# ancil_soil_albedo

The `ancil_soil_albedo` app generates the soil albedo ancillary (STASH `m01s00i216`) on the target model grid, and makes it consistent with the land cover type fraction (LCT) ancillary produced by `ancil_lct`. The precise steps performed by the task are:

1. The source soil albedo dataset is regridded onto the LCT grid using a two-stage regridding scheme.
2. Any grid cell where the LCT ice tile fraction is 0.5 or greater has its soil albedo set to a fixed value of 0.75.
3. Remaining coastline mismatches between the regridded soil albedo and the LCT land-sea mask are filled using nearest-neighbour spiral search, so that soil albedo values exist everywhere the LCT considers land.

Note that the source cannot be made consistent with soil vs. non-soil surface types in the LCT, since it isn't known which soil albedo values correspond to which surface type — only the ice adjustment and coastline fill are applied.

## Env Arguments

#### Inputs
* `ALBEDO_SOURCE`: Source soil albedo dataset, in a format loadable by ANTS (e.g. NetCDF), providing the raw soil albedo climatology to be regridded.
* `LCT_SOURCE`: Land cover type fraction ancillary (`qrparm.veg.frac`), produced by `ancil_lct`. Provides the target grid, the surface-type `pseudo_level` coordinate used to locate the ice tile, and the ice fraction used for the albedo consistency adjustment.

#### Outputs
* `OUTPUT`: Path to write the resulting soil albedo field to, in UM PP format, plus NetCDF format with the `.nc` suffix added.

#### Parameters
* `MODEL_ICE_TILE`: Which model surface (pseudo-level) index in `LCT_SOURCE` corresponds to the ice tile, used to identify ice fraction when applying the ice albedo adjustment. Not given a default in this app's `rose-app.conf`, so must be supplied by the wider suite (as with the same-named parameter in `ancil_lct`).
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
