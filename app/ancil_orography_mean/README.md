# ancil_orography_mean

The `ancil_orography_mean` app derives the mean orography on the target grid from the (unfiltered) source orography, smoothing it as needed so that it does not exceed a maximum slope. It is used to produce a stable mean orography field ahead of the sub-grid orographic wavedrag calculation performed by `ancil_orographic_wavedrag`.

The mean orography is derived as follows:

1. The unfiltered source orographic height field is conservatively regridded onto the target grid.
2. Since orography blending is enabled, the regridded mean is filtered with the Raymond filter (see `ants.analysis.filters.raymond`) using the smallest of the configured epsilon values, giving a "fine" orography.
3. If the fine orography's maximum gradient exceeds `MAX_SLOPE`, the app searches through the remaining (increasing) epsilon values for the smallest one that produces a "smooth" orography whose maximum gradient satisfies `MAX_SLOPE`. If none of the epsilon values satisfy `MAX_SLOPE`, the largest epsilon value is used instead.
4. The fine and smooth orographies are then blended together: at each point where the fine orography's gradient exceeds `MAX_SLOPE`, the field is progressively low-pass filtered and blended towards the smooth orography, repeating for up to `MAX_SMOOTH` iterations or until no gradients remain steeper than `MAX_SLOPE`.

Blending requires at least two epsilon values to be configured (a hard-coded default of eight is used here); this is not exposed as an env var. The Raymond filter epsilon values and the choice of the maximum-slope-from-smooth strategy are also not exposed as env vars for this task and use the script's defaults.

## Env Arguments

#### Inputs
* `OROGRAPHY_SOURCE`: Unfiltered source orography dataset.
* `TARGET_LSM`: Target land-sea mask defining the destination grid (`land_binary_mask`), from `ancil_lct`. Set globally for the suite.

#### Outputs
* `OUTPUT`: Path to write the (blended, smoothed) mean orography to, in UM ancil (PP) format and NetCDF format. This intermediate mean orography (STASH `m01s00i033`) is consumed as `MEAN_OROGRAPHY_INPUT` by `ancil_orographic_wavedrag`.

#### Parameters
* `MAX_SLOPE`: Maximum permitted orography gradient. Used both to decide which smoothing epsilon to select and as the threshold for identifying regions requiring blending towards the smoothed orography.
* `MAX_SMOOTH`: Maximum number of blending/smoothing iterations to perform before giving up even if steep slopes remain.
* `ANTS_CONFIG`: ANTS config file. Typically the task's `rose-app.conf`.
