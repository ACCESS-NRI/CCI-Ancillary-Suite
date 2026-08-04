# ancil_C4_fraction

The `ancil_C4_fraction` app splits the combined C3/C4 grass fraction in a land cover type fraction ancillary into separate C3 and C4 grass classes, using the ISLSCP II C4 Vegetation Percentage (Still et al.) dataset to determine the split. The injection is performed using the relation:

```
vegfrac_C4_grass = min(C4_still, vegfrac_C3_grass)
vegfrac_C3_grass = vegfrac_C3_grass - vegfrac_C4_grass
```

where the input C3 grass fraction is assumed to represent the combined C3 and C4 grass fraction prior to the split. Before injecting the C4 data, it is regridded onto the land cover type fraction grid and any points inconsistent with the CCI land-sea mask are filled using a nearest-point fill. The app raises an error if the land cover type fraction dataset already has a non-zero C4 grass fraction, to avoid double-counting.

## Env Arguments

#### Inputs
* `VEGETATION_SOURCE`: Land cover type fraction dataset produced by `ancil_lct`, prior to the C3/C4 split. Must be a NetCDF file containing a cube with STASH code `m01s00i216`.
* `C4_SOURCE`: Regridded ISLSCP II C4 Vegetation Percentage source data, expressed as a percentage (converted to a fraction internally by dividing by 100).

#### Outputs
* `OUTPUT`: Path to write the output to, in UM PP format and also in NetCDF format.

#### Parameters
* `MODEL_C3_GRASS_TILE`: Index of the C3 grass class in the land cover type fraction dataset. Set in `flow.cylc` based on the specified land model.
* `MODEL_C4_GRASS_TILE`: Index of the C4 grass class in the land cover type fraction dataset, into which the split C4 fraction is written. Set in `flow.cylc` based on the specified land model.
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
