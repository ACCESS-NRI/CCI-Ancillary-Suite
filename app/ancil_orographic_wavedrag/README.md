# ancil_orographic_wavedrag

The `ancil_orographic_wavedrag` app calculates the sub-grid orographic (wavedrag) ancillaries needed by the UM's orographic gravity wave drag scheme: the filtered mean orography, its standard deviation within each target grid box, the sub-grid gradient correlation terms, and the derived LM97 (Lott and Miller 1997) anisotropy, orientation and slope parameters.

The methodology is:

1. The Raymond-filtered source orography (produced by `ancil_orography_filtered_preproc`) has the mean-orography gradient field removed, by linearly regridding the mean orography (from `ancil_orography_mean`) back onto the source grid and subtracting it from the filtered source, giving the "filtered source minus mean gradient" (SMG) field.
2. From the SMG, the x and y gradient components are calculated and combined into the sub-grid correlation terms `sigma_xx`, `sigma_yy` and `sigma_xy`, which are regridded onto the target grid.
3. The standard deviation of the source orography within each target grid box is calculated from the SMG and its regridded (target-grid) counterpart.
4. All output fields are set to zero over open ocean, identified via a floodfill of the target land-sea mask, and any resulting `nan` values (occurring at the poles) are set to zero.
5. The LM97 parameters — orographic sub-grid anisotropy, orientation and slope — are derived from `sigma_xx`, `sigma_yy` and `sigma_xy` (see UMDP 022, section 3.1).

This produces the filtered mean (`m01s00i033`), standard deviation (`m01s00i034`), `sigma_xx` (`m01s00i035`), `sigma_yy` (`m01s00i036`), `sigma_xy` (`m01s00i037`), and the LM97 anisotropy (`m01s06i248`), orientation (`m01s06i249`) and slope (`m01s06i250`) fields. The mean orography (after ocean masking) is separately re-saved to its own output path. The task requires that its two output paths differ, and will raise an error otherwise.

## Env Arguments

#### Inputs
* `OROGRAPHY_SOURCE`: Raymond-filtered source orography, from `ancil_orography_filtered_preproc`. Must contain a `surface_altitude_filtered` field.
* `MEAN_OROGRAPHY_INPUT`: Mean orography on the target grid, from `ancil_orography_mean`. Must contain a field with STASH `m01s00i033`.
* `TARGET_LSM`: Target land-sea mask defining the destination grid (`land_binary_mask`), from `ancil_lct`. Set globally for the suite.

#### Outputs
* `MEAN_OROGRAPHY_OUTPUT`: Path to write the ocean-masked mean orography to, in UM ancil (PP) format and NetCDF format.
* `TARGET_OROGRAPHIC_LSM`: Path to write the standard deviation, sigma and LM97 wavedrag fields to, in UM ancil (PP) format and NetCDF format. Set globally in `flow.cylc` rather than this app's own `rose-app.conf`, and consumed under this same name by downstream tasks such as `regrid_aerosol` and `regrid_ozone`.

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
