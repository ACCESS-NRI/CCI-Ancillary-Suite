# CCI-Ancillary-Suite

## About

A [Cylc8](https://cylc.github.io/cylc-doc/stable/html/index.html) workflow, utilising the [Rose plugin](https://metomi.github.io/rose/doc/html/index.html), for generating ancillaries for land and atmosphere land models, targeting ACCESS3 models e.g. CABLE, ACCESS-AM3. Based on the 300m resolution CCI Land Cover dataset.

## Usage

Configure the workflow by modifying the `rose-suite.conf` file. The top level configuration options are described in the [Configuration][#top-level-configuration] section below. The workflow by default runs and uses storage under the user's default project `$PROJECT`. The working directory is `/scratch/<project>/<user>/cylc-run/CCI-Ancillary-Suite/<runID>/, with the output data located in that directory under `share/data`.

## Requirements

To run the suite on Gadi, membership is required to the following projects:

* hr22
* vk83
* access
* cm45

## Top Level Configuration

The top level configuration options, located in `rose-suite.conf`, are:

* `COMPUTE_PROJECT`: Which project to request compute resources from. Defaults to `${PROJECT}`.
* `STORAGE_PROJECT`: Which project to use for storage. Becomes the `<project>` in the above working directory path. Defaults to `${PROJECT}`.
* `LAND_MODEL`: Which land model to target. Affects the generation of vegetation ancillaries, which depend tiles active in the model. Possible options are `CABLE` and `JULES`. Defaults to `CABLE`.
* `GRID_SOURCE`: How the land/sea mask should be defined for the workflow. There are 3 possible options:
    * `ocean_mesh`: Generate a land/sea mask from a provided ocean mesh, then project the land cover onto it. If this is supplied, the `RESOLUTION` option must also be defined.
    * `land_cover`: Use the provided land cover dataset to determine the land/sea mask.
    * `land_sea_mask`: Use a supplied land.sea mask and project the land cover onto it.
    Defaults to `ocean_mesh`.
* `GRID_FILE`: The file defining the grid to use. The form of this file depends on the `GRID_SOURCE`. When the `GRID_SOURCE` is:
    * `ocean_mesh`, the `GRID_FILE` should be a CM3 ocean mesh file.
    * `land_cover`, the `GRID_FILE` should be a UM grid namelist file.
    * `land_sea_mask`, the `GRID_FILE` should be a land sea mask in UM ancillary format.
    Defaults to `"/scratch/rp23/lw5085/access-om3-025deg-ESMFmesh.nc"`.
* `VERTICAL_DISCRETIZATION`: Vertical discretization definition in UM namelist format. Defaults to an 85 level discretization at `"/g/data/access/TIDS/UM/ancil/data/namelist/vertical/vert_85/latest/vert_85"`.
* `CALENDAR`: Which calendar to use for the time-dependent ancillaries. Possible options are `360day` and `gregorian`. Defaults to `360day`.
* `RESOLUTION`: Used only when `GRID_SOURCE="ocean_mesh"`. Must be a valid UM resolution specifier e.g. `n96e`, `n512e`.

## Common App Configurations

There are some apps that are likely to be modified for scientific exploration:

* `ancil_lct`: The application which maps the CCI land cover classes to land cover classes for the target land model. Modify the `TRANSFORM_FILE` environment variable in the app configuration to use a different mapping from CCI land cover classes to the model's classes. An absolute path can be supplied, or the new transform can be added to the app's `file` directory.
* `ancil_LAI`: The application which computes the leaf area index (LAI) per plant functional type by combination of a source LAI dataset and the generated land cover map for the target model. The weightings applied to distribute the LAI between classes is defined by a weightings file, specified by the `RELATIVE_WEIGHTS` environment variable in the app configuration.
* `ancil_canopy_height`: The application which computes the canopy height per plant functional type by combination of a source canopy heights dataset and the generated LAI. A set of height factors determine how canopy heights are computed, defined the file specified by the `HEIGHT_FACTORS` environment variable in the app configuration.
