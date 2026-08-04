# ancil_split_grass_types

The `ancil_split_grass_types` is a task intended to overcome a shortcoming of the original CCI land cover dataset, which is its limited classification of grassy/short vegetation. An NCAR dataset, which splits grasses into C3, C4 and arctic grasses as well as different shrub types is used to determine how the CCI classifications should be split.

The short vegetation types in the input surface area fractions are grouped together, and then re-distributed by the NCAR dataset. Effectively, the original land cover is used to determine how much short vegetation there is, and the NCAR dataset determines how that short vegetation should be split into the various sub-types of short vegetation. In the CABLE example, the short vegetation types are shrub, C3 grass, C4 grass and tundra.

This relies on a non-provenanced dataset, and should be retired when an appropriate replacement dataset that is published is found.

## Env Arguments

#### Inputs
* `GRASS_SOURCE`: NCAR dataset containing more detailed short vegetation classifications.
* `ORIG_FRACTIONS`: Original surface fractions to be processed.
* `VEGETATION_MAPPING`: JSON file describing the mapping from `GRASS_SOURCE` vegetation types to `ORIG_FRACTIONS` vegetation types.

#### Outputs
* `OUTPUT`: Path to write the output to in PP format and NetCDF format with the `.nc` extension.
