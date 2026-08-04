# ancil_aeroclim_preproc

The `ancil_aeroclim_preproc` app pre-processes a set of GLOMAP MODE aerosol source PP files, produced by a previous Unified Model run, into a monthly climatology suitable for use as a GLOMAP MODE ancillary. For each of the 12 months, it loads all matching PP files, computes the time mean of each field, and updates the resulting cube metadata so that it is compatible with what GLOMAP MODE expects in the UM and LFRic:

1. Load all PP files matching a given month from the source directory.
2. Collapse each field over the time coordinate to produce a monthly mean, realising the data during this step for performance.
3. Rewrite the STASH attribute, derive the required `stashcode` attribute, and rename the NetCDF variable to the name expected by the LFRic interface.
4. Convert the time coordinate to hours since 1970-01-01, using either a Gregorian or 360-day calendar, and set monotonic time bounds spanning the full climatology.
5. Merge all months/fields together and write the climatology out using the ANTS UKCA NetCDF saver.

## Env Arguments

#### Inputs
* `AEROCLIM_SOURCE`: Directory/path prefix used to build a glob pattern (`{AEROCLIM_SOURCE}/*{month}.pp`) that locates the source PP files for each month.

#### Outputs
* `OUTPUT`: Path to write the pre-processed aerosol climatology to, in NetCDF format (written using the ANTS UKCA saver, `ants.io.save.ukca_netcdf`).

#### Parameters
* `METADATA_SOURCE`: Free-text description of the provenance of the source data (e.g. which UM experiment and time period it was generated from). Written into the `source` attribute of every output cube.
* `CALENDAR`: Calendar to use for the output time coordinate, either `gregorian` or `360day`. Set in `rose-suite.conf`.
