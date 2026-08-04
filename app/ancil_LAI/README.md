# ancil_LAI

The `ancil_LAI` app derives the leaf area index (LAI) for each plant functional type (PFT) of the land cover type fraction ancillary, from a total (all-vegetation) LAI source. The total LAI is related to the per-PFT LAI through `L_tot = sum_i(L_i * f_i)`, where `f_i` is the fraction of land covered by PFT `i`. Assuming the LAI of each PFT is fixed relative to the others (via a supplied set of relative weights), this system can be solved by substitution to give each PFT's LAI at every grid point. The precise steps performed by the task are:

1. Regrid the total LAI source onto the grid of the land cover type fraction ancillary.
2. Where the total LAI is missing across some (but not all) months, set it to a floor value of `1e-10`; any other unmasked value below this floor is also raised to it.
3. At each grid point, solve for the per-PFT LAI via substitution, using the supplied relative weights, wherever the total vegetation fraction across the functional types is non-zero.
4. Where the total vegetation fraction of the functional types is non-zero but less than 20% of the grid box, mask the resulting LAI PFTs, to avoid unreasonably high LAI values from being resolved.
5. Fill the masked/missing points using a nearest neighbour (spiral) search, so that the result is consistent with the coastlines of the land cover type fraction ancillary.

## Env Arguments

#### Inputs
* `LAI_SOURCE`: Total (all-PFT) LAI on the target grid. Must be a NetCDF dataset containing a `leaf_area_index` variable.
* `LAND_COVER`: Land cover type fraction ancillary. Must be a NetCDF dataset with a cube of STASH code `m01s00i216`, used to define the target grid, the functional type fractions and the land-sea mask.

#### Outputs
* `OUTPUT`: Path to write the output to, in UM PP format and also in NetCDF format.

#### Parameters
* `RELATIVE_WEIGHTS`: JSON file describing the relative LAI weighting between PFTs (`relative_weights` and `jules_classes` keys, one weight per class). Defaults to a CABLE or JULES specific set of weights depending on the chosen land model; can be replaced with a custom set of weights.
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
