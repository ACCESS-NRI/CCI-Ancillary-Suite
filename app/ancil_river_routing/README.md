# ancil_river_routing

The `ancil_river_routing` app derives the river routing sequence and direction ancillary (`qrparm.rivseq`), which describes how runoff is routed across each land grid box towards the ocean. The direction field encodes, for each point, which of its eight neighbours (or an ocean/inland outflow) the flow proceeds to, using values 1-8 for compass directions, 9 for an outflow pour point into the ocean, and 10 for an inland basin. The precise steps performed by the task are:

1. The land cover fraction (from `ancil_lct`) is regridded onto the river routing grid using an area-weighted mean.
2. The regridded land cover fraction is used to identify points that are entirely ocean, and coastal points that have no river flow going into them.
3. At points that are entirely ocean, and coastal points with no inflow, the river direction and sequence are set to missing data.
4. At coastal points that do have flow going into them, the river direction is set to 9 (a pour point into the ocean).
5. Optionally (via `--make_nemo_rivers`), rivers are sorted into a unique order and UM/NEMO river number ancillary files are produced using the ocean runoff and ORCA domain data. This flag is not passed by the app's command, so this step is not currently exercised in this suite.
6. The output grid definition is inherited from the land cover fraction input, since the river routing source grid is otherwise decoupled from the UM grid.

The land threshold below which a point is treated as an outflow point is fixed at 0.5 in the app's command, appropriate for GC5-and-later style ancillaries (allowing outflow points to fall on the coast, rather than requiring a point to be entirely sea as with GC4 and earlier). A comment in the app's `rose-app.conf` notes uncertainty over whether the NEMO river routing outputs are actually needed for AM3.

## Env Arguments

#### Inputs
* `ROUTING_SOURCE`: Source river routing sequence and direction data, in NetCDF format, produced by `ancil_river_routing_preproc`.
* `LANDFRAC`: Land area fraction ancillary (`land_area_fraction`), produced by `ancil_lct`, used to determine ocean/coastal points on the river routing grid.
* `RUNOFF_SOURCE`: Ocean runoff dataset used for ocean-only runs, only consumed when NEMO river number generation (`--make_nemo_rivers`) is enabled, which it currently is not in this app's command.
* `ORCA_DOMAIN`: NEMO domain file containing the land-sea mask for the NEMO model on the ORCA grid, likewise only consumed when NEMO river number generation is enabled.

#### Outputs
* `OUTPUT`: Path to write the river routing sequence and direction fields to, in UM PP format, plus NetCDF format with the `.nc` suffix added.

#### Parameters
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
