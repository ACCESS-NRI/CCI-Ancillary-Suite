# ancil_soil_roughness

The `ancil_soil_roughness` app derives the soil (bare ground) roughness length ancillary from a source roughness climatology, filtered and blended using a leaf area index (LAI) climatology. The precise steps performed by the task are:

1. Compute the time-mean and time-minimum of the supplied LAI climatology, and regrid these to the source roughness grid.
2. Mask out source roughness values that are considered unreliable due to vegetation cover: points where the roughness is greater than 0.08 m and the time-mean LAI is greater than 4, or where the time-minimum LAI is greater than 2.
3. Regrid the (masked) soil roughness source to the target grid, using a bespoke area-weighted regrid of the roughness converted to a flux-weighted quantity assuming a blending height of 10 m.
4. Combine the raw, masked and gap-filled regridded fields: ocean points are masked, land points with no data are set to a fill value (1.0e-3), and points with sufficient "known good" (non-masked, non-high-roughness) source data take precedence over the fill value.
5. Any points still missing after this process (i.e. missing from the raw source data itself, not simply filtered by the LAI mask) are filled in using a spiral search interpolation, consistent with the target land-sea mask.

Note that the blending height of 10 m is currently hardcoded rather than derived from the source roughness dataset.

## Env Arguments

#### Inputs
* `ROUGHNESS_SOURCE`: Source roughness length dataset. Must be readable as a cube with `STASH` `m01s00i097`.
* `LAI_SOURCE`: Leaf area index source. Should be the intermediate pre-processed LAI climatology file (i.e. not yet regridded to the target grid, and not split onto land surface tiles).
* `TARGET_LSM`: Target land-sea mask defining the domain and grid to output the soil roughness on.

#### Outputs
* `OUTPUT`: Path to write the output to, in NetCDF format and, unless run with `netcdf_only`, also in UM PP/ancillary format.

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
