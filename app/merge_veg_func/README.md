# merge_veg_func

The `merge_veg_func` app merges the generated leaf area index (LAI) and canopy heights datasets into a single vegetation function ancillary. It uses ANTS's generic `ancil_2anc` file-format translation utility, which loads one or more cubes from the supplied source file(s) and writes them to the output ancillary/NetCDF, assuming all fields present in the input are to be included in the output and that their metadata is already complete and accurate (the `history` attribute is stripped before saving). No other processing or regridding is performed.

## Env Arguments

#### Inputs
* `LAI_SOURCE`: Generated leaf area index dataset to merge.
* `CANOPY_HEIGHTS_SOURCE`: Generated canopy heights dataset to merge.

#### Outputs
* `OUTPUT`: Path to write the merged vegetation function ancillary to, in NetCDF format and UM PP/ancillary format.
