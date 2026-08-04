# ancil_lct

The `ancil_lct` app converts a source land cover dataset, which characterises areas of land into distinct land cover classes, to a set of surface area fractions usable by land models. It also applies some post-processing steps to the ice fractions to be in line with what JULES expects from the land fractions. The precise steps performed by the task are:

1. Convert the land cover classes to land area fractions on the specified grid using a supplied mapping file, which describes how the original land cover classes correspond to the target model's cover classes.
2. Any grid cells fractions with ice fraction greater than 0.5 are set to 1.0, and any less are set to 0.0 and then the fractions are re-normalised to 1.0. This is required as JULES has special handling for ice "soil" columns, and shared columns between tiles.
3. Non-glacial ice is removed. Non-glacial ice is defined by ice points that are "isolated". Any non-glacial ice is replaced with bare soil.

## Env Arguments

#### Inputs
* `LAND_COVER_SOURCE`: Source land cover dataset. Must be a NetCDF dataset with a `land_cover_lccs` variable which has `flag_values` and `flag_meanings` attributes.
* `TRANSFORM_FILE`: JSON file which provides a mapping from the land cover dataset to the target model's surface types. The required entries in the JSON file are:
  - `cover_map`: A matrix (list of lists) describing how to apportion the source land cover types to the model's surface types. Should be of size `N_source_classes x (1+N_model_types + 1)`, where the `1+` is for the ocean (the ocean should be the first index, with model surface types following).
  - `source`: List of source class names. Should correspond with entries in the `flag_meanings` attribute of `land_cover_lccs`.
  - `target`: List of target model surface indices, with 0 denoting ocean.
* `GRID_FILE`: A UM grid definition namelist to use to define the domain. Cannot be included with `TARGET_LSM`.
* `TARGET_LSM`: A land-sea mask to use to defined the domain. Cannot be included with `GRID_FILE`.

#### Outputs
* `OUTPUT`: Path to write the output to. Writes to the provided path in UM PP format, and also to NetCDF format with the `.nc` suffix added.
* `LSM_OUTPUT_PATH`: When using `GRID_FILE` to define the domain, the task also writes out the land-sea mask and land fractions into this directory in PP and NetCDF format.

#### Parameters
* `MODEL_SOIL_TILE`: Which model surface index corresponds to the soil type, used when removing the non-glacial ice.
* `MODEL_ICE_TILE`: Which model surface index corresponds to the ice type, used when performing ice adjustments.
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
