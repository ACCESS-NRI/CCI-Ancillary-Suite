# merge_soils

The `merge_soils` app merges the previously generated soil hydrology and soil albedo ancillaries into a single soils ancillary file.

This is required because the UM Continuous Ancillary Processing (CAP) reads specifically the first 10 fields from the input file when generating the snow ancillaries (`ancilSmcsnow`), and there are more than 10 fields present across the combined soils ancillary once hydrology and albedo are merged, meaning `ancilSmcsnow` is dependent on field order. The task therefore:

1. Loads the hydrology and albedo source cubes.
2. Removes any scalar `time` coordinate present on the cubes, so that all cubes share consistent coordinates.
3. Removes the `history` attribute from each cube.
4. Sorts the cubes into a fixed STASH item order (`40, 41, 43, 207, 47, 44, 46, 48, 220, 223, 8, 418, 419, 420`), moving the fields that CAP does not need for `ancilSmcsnow` generation to the end of the file, to avoid a segmentation fault.
5. Saves the combined result to output.

## Env Arguments

#### Inputs
* `HYDROLOGY_SOURCE`: Generated soil hydrology ancillary (from the `ancil_soil_hydrology`/`merge_soils` preprocessing chain) to merge.
* `ALBEDO_SOURCE`: Generated soil albedo ancillary to merge.

#### Outputs
* `OUTPUT`: Path to write the merged soils ancillary to, in UM PP/ancillary format and NetCDF format.
