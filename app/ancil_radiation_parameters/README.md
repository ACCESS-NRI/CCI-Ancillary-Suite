# ancil_radiation_parameters

The `ancil_radiation_parameters` app calculates the unfiltered mean orography and its gradient fields on the target grid, for use by the UM's radiation scheme. It should be used in combination with a source orography that has not undergone Raymond filter preprocessing (i.e. the output of `ancil_orography_unfiltered_preproc`).

The methodology is:

1. The mean orography is derived by conservatively regridding the source orographic height field onto the target grid.
2. The x and y gradient (slope) fields are calculated from the source orographic height field and conservatively regridded onto the target grid. Points beyond a fixed polar cutoff of 88 degrees are set to zero, as are any `nan` points.
3. All fields are set to zero over open ocean, identified via a floodfill of the supplied target land-sea mask.

This generates the unfiltered orography (`m01s00i007`), slope_x (`m01s00i005`) and slope_y (`m01s00i006`) fields.

## Env Arguments

#### Inputs
* `OROGRAPHY_SOURCE`: Unfiltered source orography, from `ancil_orography_unfiltered_preproc`. Must contain an `unfiltered_surface_altitude` field.
* `TARGET_LSM`: Target land-sea mask defining the destination grid (`land_binary_mask`), from `ancil_lct`. Set globally for the suite.

#### Outputs
* `OUTPUT`: Path to write the unfiltered orography and slope fields to, in UM ancil (PP) format and NetCDF format.

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
