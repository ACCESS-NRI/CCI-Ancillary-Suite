# ancil_orographic_formdrag

The `ancil_orographic_formdrag` app generates the orographic formdrag ancillaries used by the UM: half of the peak-to-trough height (`m01s00i018`) and the silhouette of orography per unit area (`m01s00i017`). These fields quantify sub-grid orographic roughness and are used by the atmosphere model's orographic drag scheme. The precise steps performed by the task are:

1. Regrid the silhouette and half-height source fields (produced by `ancil_orographic_formdrag_preproc`) onto the target grid using a two-stage regridding scheme, and make the result consistent with the supplied land-sea mask.
2. Set all ocean points to zero, rather than leaving them as missing data. Leaving them as missing data causes the UM to crash on the first time step.
3. Divide the half-height field by `2*sqrt(2)` so that it is consistent with the equivalent field produced by the CAP.

## Env Arguments

#### Inputs
* `OROGRAPHY_SOURCE`: Preprocessed orography source, from `ancil_orographic_formdrag_preproc`. Must contain the silhouette (STASH `m01s00i017`) and half-height (STASH `m01s00i018`) fields.
* `TARGET_LSM`: Target land-sea mask defining the destination grid, from `ancil_lct`. Must contain a field with STASH `m01s00i030`.

#### Outputs
* `OUTPUT`: Path to write the output to, in UM ancil (PP) format and NetCDF format.

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
