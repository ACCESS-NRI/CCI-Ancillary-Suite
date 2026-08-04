# ancil_soil_hydrology

The `ancil_soil_hydrology` app derives the soil hydrology, thermal and carbon parameters required by the land model from a soils source field and a Cosby-parameter lookup table. The precise steps performed by the task are:

1. For each source location, look up the possible soil types (and their area share) associated with the `MU_GLOBAL` soil unit identifier, using the supplied JSON lookup table.
2. Determine the dominant soil type for each target grid cell, using ESMF-derived conservative area weights between the source grid and the target grid.
3. Look up the sand-silt-clay based (Clapp-Hornberger) hydrology, carbon and thermal property parameters corresponding to the dominant soil type for each target cell.
4. Make the result consistent with the ice fraction of the land cover type ancillary: locations where every parameter is zero are masked, a fill search is used to make the parameter mask consistent with the land-sea/ice mask, and parameter values at ice locations are then overridden (all parameters set to 0, except soil thermal capacity and soil thermal conductivity, which are set to fixed values of 630000 and 0.2650 respectively).
5. Save the resulting soil hydrology parameters (Brooks-Corey exponent, saturated soil suction, saturated hydraulic conductivity, and volumetric soil moisture at saturation/wilting/critical points), bulk density, soil carbon content and thermal properties (conductivity and capacity). Volume fractions of sand, silt and clay are also saved for reference, though these are not used directly by the model.

## Env Arguments

#### Inputs
* `SOILS_SOURCE`: Source soils dataset, providing a geo-located `MU_GLOBAL` soil unit field.
* `LCT_SOURCE`: Land cover type fraction ancillary (from the `ancil_lct` app), used to identify ice-covered locations for consistency checking. Must be readable as a cube with `STASH` `m01s00i216`.
* `SOILS_LOOKUP`: JSON soil lookup table mapping `MU_GLOBAL` values to unique soil IDs (with their area share), and mapping each unique soil ID to its sand-silt-clay based Cosby hydrology/thermal/carbon parameters.

#### Outputs
* `OUTPUT`: Path to write the output to. The app always writes both UM ancillary (PP) format and NetCDF format (with the `.nc` suffix added), since the `--use-new-saver` flag is unconditionally passed in the command.

#### Parameters
* `MODEL_ICE_TILE`: Which model surface index corresponds to the ice tile in `LCT_SOURCE`, used to identify ice-covered locations so the derived soil parameters can be overridden consistently with the ice mask.
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
