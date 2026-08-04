# ancil_orography_filtered_preproc

The `ancil_orography_filtered_preproc` app produces a Raymond-filtered version of the source orography, for use as an input to `ancil_orographic_wavedrag`. An isotropic Raymond filter with a filter length scale of 6 km is applied to the source surface altitude field to remove sub-6 km scale features, any resulting negative heights are clipped to zero, and the field is renamed to `surface_altitude_filtered`.

This preprocessing step does not depend on the model resolution, target grid or any other ancillary, so it can in principle be run independently of the rest of the suite.

## Env Arguments

#### Inputs
* `SOURCE`: Source orography dataset to be filtered.

#### Outputs
* `OUTPUT`: Path to write the filtered `surface_altitude_filtered` field to, in NetCDF format only.
