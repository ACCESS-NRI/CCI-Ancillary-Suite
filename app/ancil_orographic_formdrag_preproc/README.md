# ancil_orographic_formdrag_preproc

The `ancil_orographic_formdrag_preproc` app corrects known metadata problems in the source GLOBE ash file used to derive the orographic formdrag ancillaries, ahead of processing by `ancil_orographic_formdrag`. The source file has incorrect values set for grid staggering, `lbcode` and the missing data indicator. The following corrections are applied to every field:

1. `grid_staggering` is set to 6.
2. `bmdi` is set to the mule real missing-data indicator.
3. `lbcode` is set to 1.
4. The silhouette field (STASH `m01s00i017`) is multiplied by a factor of 1.823310 to account for the reduced resolution of the globe30 source data. This scaling is only appropriate for globe30 source data. The half-height field (STASH `m01s00i018`) is passed through unchanged.

The app can only process a single source file per invocation, and will raise an error if a field is encountered with a STASH code other than 17 or 18.

## Env Arguments

#### Inputs
* `OROGRAPHY_SOURCE`: Source orographic ash file (UM ancil format) containing the silhouette (STASH `m01s00i017`) and half-height (STASH `m01s00i018`) fields to be corrected.

#### Outputs
* `OUTPUT`: Path to write the corrected fields to, in UM ancil format.

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
