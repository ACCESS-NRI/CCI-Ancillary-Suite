# preprocess_CCI_water

The `preprocess_CCI_water` app merges the CCI land cover dataset for a given year (300 m resolution) with the CCI permanent water bodies dataset (150 m resolution), so that the water bodies dataset can be treated as ground truth for identifying ocean and inland water, overriding the (coarser and less reliable) water classification in the land cover dataset. This is an optional pre-processing step that should only be run when a new year of the CCI land cover dataset needs to be processed. The task:

1. Combines and upscales the two source water bodies datasets (an ocean-only mask and an all-water mask) from 150 m to 300 m resolution, using a windowed majority-vote scheme, to produce a single water classification: ocean, land, or inland water.
2. Overlays this water classification onto the land cover data: the original generic `water` land cover class is replaced with two new classes, `sea_ocean_water` and `inland_water`, assigned according to the upscaled water bodies dataset. Any cells still flagged as generic water by the land cover dataset, but not confirmed as water by the permanent water bodies dataset, are filled in with the nearest valid land classification (nearest point considered land by the permanent water dataset).
3. Drops land cover ancillary variables that are not needed downstream, to reduce disk usage.
4. Adds a `latitude_longitude` coordinate system, as required by ANTS/Iris to correctly interpret the data.
5. Corrects assorted metadata issues in the source land cover dataset (an unused `no_data` flag, `valid_min`/`valid_max` set for unsigned data despite the flags being signed, and an incorrect `_Unsigned` encoding attribute).
6. Writes the corrected, merged dataset to NetCDF.

Note that the water bodies source files are expected to have fixed names (`ESACCI-LC-L4-WB-Ocean-Map-150m-P13Y-2000-v4.0.tif` and `ESACCI-LC-L4-WB-Map-150m-P13Y-2000-v4.0.tif`) within `WATER_BODIES_SOURCE`.

## Env Arguments

#### Inputs
* `LAND_COVER_SOURCE`: Source CCI land cover dataset for the target year (`CCI_YEAR`). NetCDF dataset with an `lccs_class` variable carrying `flag_values`/`flag_meanings` attributes.
* `WATER_BODIES_SOURCE`: Directory containing the CCI permanent water bodies source datasets (fixed-name GeoTIFFs for the ocean-only and all-water masks).

#### Outputs
* `OUTPUT`: Path to write the merged land cover and water bodies dataset to, in NetCDF format.

#### Parameters
* `CCI_YEAR`: Which year of the CCI land cover dataset to target; used to construct the default `LAND_COVER_SOURCE` path.
