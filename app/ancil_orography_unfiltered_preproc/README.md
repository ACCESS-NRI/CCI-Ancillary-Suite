# ancil_orography_unfiltered_preproc

The `ancil_orography_unfiltered_preproc` app converts the source orography dataset to NetCDF, without applying any filtering, for use as an input to `ancil_radiation_parameters` and (via `ancil_orography_mean`) `ancil_orographic_wavedrag`. The field is renamed to `unfiltered_surface_altitude` and given STASH code `m01s00i007`. Converting to NetCDF at this stage allows the data to be realised in sub-sets by the downstream processing tasks.

This preprocessing step does not depend on the model resolution, target grid or land cover, so it can in principle be run independently of the rest of the suite.

## Env Arguments

#### Inputs
* `SOURCE`: Source orography dataset.

#### Outputs
* `OUTPUT`: Path to write the `unfiltered_surface_altitude` field to, in NetCDF format only.
