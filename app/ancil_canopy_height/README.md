# ancil_canopy_height

The `ancil_canopy_height` app derives canopy height ancillaries from the leaf area index (LAI), then overrides the tree plant functional types (PFTs) using an independent tree height dataset. The precise steps performed by the task are:

1. Calculate a first-pass canopy height for every PFT from the LAI using the relation `canopy_height = height_factor * LAI^(2/3)`, where the height factor is looked up per-PFT from a supplied JSON mapping.
2. Regrid the Simard/Pinto global vegetation tree height dataset onto the target grid.
3. Make the regridded tree field consistent with the land-sea mask of the LAI by performing an index-based nearest neighbour search, constraining the search in `y` to a distance of 500 km.
4. Override the tree PFTs calculated in step 1 with the values derived from the tree dataset.
5. Since the canopy heights are time-invariant but must be stored alongside the periodic 12-month LAI fields, the final field is the maximum canopy height across the year, duplicated across all 12 months.

Note that this app is run with the `--netcdf-only` flag, so the output ancillary is only written in NetCDF format (no UM PP/ancil file is produced).

## Env Arguments

#### Inputs
* `LAI_SOURCE`: Leaf area index ancillary (NetCDF, cube with STASH code `m01s00i217`) used to derive the initial canopy heights.
* `TREE_HEIGHT_SOURCE`: Simard/Pinto 3D Global Vegetation canopy height dataset, used to override the tree PFTs.

#### Outputs
* `OUTPUT`: Path to write the canopy height ancillary to, in NetCDF format only.

#### Parameters
* `HEIGHT_FACTORS`: JSON file mapping a canopy height factor to each JULES/model class (`canopy_height_factors` and `jules_classes` keys). Defaults to a CABLE or JULES specific set of factors depending on the chosen land model; can be replaced with a custom set of factors.
* `MODEL_TREE_TILES`: Comma-separated list of PFT indices considered "tree" types, which are overridden using the tree height dataset. Set in `flow.cylc` based on the specified land model.
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
