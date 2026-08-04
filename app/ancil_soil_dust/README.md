# ancil_soil_dust

The `ancil_soil_dust` app generates the soil dust parent ancillary, representing the relative mass fraction of soil particles (`M_rel`) in each of 6 UKMO particle-size divisions, as required by the dust emissions scheme. It is derived from an observation-based soil dataset that only provides the clay, silt and sand mass fractions of the soil.

To relate these three broad fractions to the 6 model size divisions, the app assumes a specific form for the particle-size distribution in `dM/dlog(R)` vs `log(R)` space (following Woodward, 2001, https://doi.org/10.1029/2000JD900795): a fixed clay/silt boundary at 1 micron radius and silt/sand boundary at 25 microns, a uniform ("histogram-like") distribution across the silt and sand ranges, and a linear `dM/dlog(R)` distribution for clay that is continuous with the silt distribution at the clay/silt boundary. Given the clay and silt fractions, this defines the full distribution, which is then numerically integrated across each of the 6 model size bins (further subdivided into 10 sub-bins each for the integration) to obtain the mass fraction represented by each bin. The calculation is skipped at ocean points and at points identified as land ice.

The resulting 6 division fields are written out alongside copies of the input sand, silt and clay fraction fields.

## Env Arguments

#### Inputs
* `SOIL_SOURCE`: Soil ancillary (`qrparm.soil.nc`) containing the sand (STASH `m01s00i420`), silt (STASH `m01s00i419`) and clay (STASH `m01s00i418`) mass fraction fields.
* `LCT_SOURCE`: Land cover type fraction ancillary (`qrparm.veg.frac`), produced by `ancil_lct`, used to derive the land-sea mask and the land-ice mask. The land-ice mask is taken from the last surface type (pseudo-level) in the ancillary, which is assumed to be the ice tile — unlike `ancil_soil_albedo`, the ice tile index is not exposed as a configurable parameter here.

#### Outputs
* `OUTPUT`: Path to write the 6 soil dust mass fraction division fields (plus the input sand, silt and clay fractions) to, in UM PP format, plus NetCDF format with the `.nc` suffix added.

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
