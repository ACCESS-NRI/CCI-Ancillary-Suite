# ancil_topographic_index

The `ancil_topographic_index` app derives the mean and standard deviation topographic index ancillaries (used by JULES' TOPMODEL-based hydrology) on the target grid, and ensures they are consistent with the land cover type fraction ancillary's land-sea mask and ice fraction. The precise steps performed by the task are:

1. Area-weighted regrid the source topographic index data onto the grid defined by the land cover type fraction ancillary, to produce the mean field.
2. Derive the standard deviation field from the source data and the mean field just derived.
3. Fill all non-ice locations that have missing/zero topographic index values using a nearest-point fill, so that both fields are consistent with the land-sea mask.
4. Explicitly set the topographic index (mean and standard deviation) to 0 at ice locations, rather than leaving them masked.

Note that GA7-generation ancillaries used a different (Moore neighbourhood based) algorithm to make the fields consistent with the land-sea mask; this app uses ANTS' `FillMissingPoints` instead.

## Env Arguments

#### Inputs
* `TOPOGRAPHY_SOURCE`: Source topographic index dataset to regrid onto the target grid.
* `LCT_SOURCE`: Land cover type fraction ancillary produced by `ancil_lct`. Must be a NetCDF dataset with a cube of STASH code `m01s00i216`, used to define the target grid, land-sea mask and ice fraction for consistency handling.

#### Outputs
* `OUTPUT`: Path to write the output to, in UM PP format and also in NetCDF format. Contains both the mean (STASH `m01s00i274`) and standard deviation (STASH `m01s00i275`) topographic index fields.

#### Parameters
* `MODEL_ICE_TILE`: Index of the ice class in the land cover type fraction ancillary, used to identify ice locations that should be excluded from the fill and forced to 0.
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
