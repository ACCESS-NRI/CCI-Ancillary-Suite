# landseamask_from_ocean_mesh

The `landseamask_from_ocean_mesh` app derives the atmosphere land-sea mask and land area fractions from an ACCESS-OM3 ocean model mesh, rather than from a land cover dataset. This is used when the atmosphere land-sea mask needs to be made consistent with a coupled ocean model's grid. The precise steps performed by the task are:

1. Load the ocean mesh (an ESMF unstructured mesh file) and its element ocean/land mask.
2. Construct an equivalent unstructured UM atmosphere mesh for the requested resolution (parsed from a string such as `n96e`), including element connectivity, coordinates and areas.
3. Conservatively regrid the ocean mask from the ocean mesh onto the atmosphere mesh elements using ESMF, giving an ocean fraction per atmosphere cell, from which the land fraction is computed as `1 - ocean_fraction`. Land fractions below 0.01 are clipped to 0, and any spurious values above 1.0 are clipped to 1.0.
4. Save the land area fraction (`qrparm.landfrac`), the land binary mask (`qrparm.mask`), and a separate sea binary mask (`qrparm.mask_sea`, needed because the land mask alone would mask out fractional ocean points and break downstream ocean ancillary generation), each in both NetCDF and UM ancillary (PP) format.

Note that the `--atm-resolution` string must match the `n<number>e` naming convention used for UM ENDGame grids (e.g. `n96e`, `n512e`), and the resolution number must be even.

## Env Arguments

#### Inputs
* `OCEAN_MESH`: ACCESS-OM3 ocean mesh file to derive the land-sea mask and land fractions from. Must be an ESMF mesh-format NetCDF file with an `elementMask` variable.

#### Outputs
* `OUTPUT`: Directory to write the land-sea mask and land fraction outputs to (`qrparm.landfrac`, `qrparm.mask`, `qrparm.mask_sea`), in both NetCDF and UM ancillary (PP) format.

#### Parameters
* `RESOLUTION`: String describing the target atmosphere resolution (e.g. `n96e`, `n512e`), used to construct the atmosphere grid/mesh that the ocean mesh is regridded onto.
