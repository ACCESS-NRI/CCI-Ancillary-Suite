#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
Land cover type fraction application
************************************

Derive the land cover type fraction, land cover fraction and optionally the
landsea mask ancillaries by utilising the provided landcover types source with
its corresponding transform (crosswalk table), mapping the source classes to
the output JULES target classes.
If providing a pre-existing landsea mask, the land cover type fraction will be
made consistent with it, otherwise, the land cover type fraction field itself
will be used to derive and output a landsea mask.

Here are the steps taken:

- Transform source land classifications to the target classifications
  (ordinarily this target classification will be JULES classes).
- Regrid these target classes to the provided target grid.
- Extract the land cover fraction field ('m01s00i505').
- Set ice to whole fractions by applying a >= 50% threshold of ice to 1,
  otherwise set to 0 and re-distributing differences amongst non-zero
  land fractions (the UM defines the requirement for whole fraction ice).
- Remove non-glacial ice (driven by UM behaviour with isolated ice points).
    - Identify ice seed points by those locations of ice which are surrounded
      by 5x5 cells of ice.
    - Find contiguous ice regions using these seed points by search
      neighbouring cells (including diagonals).
    - Replace ice with soils for those locations not amongst these contiguous
      regions of ice found.
- Ensure land fractions add up to 1 and ocean fraction add up to 1 by:
    - If an optionally supplied landsea mask (``--target-lsm``) is provided, ensure
      fields are consistent with the supplied mask.
    - If the mask is to be derived by the land cover type fraction field
      (``--landseamask-out``), a threshold of >= 50% ocean is to mean ocean and
      differences are redistributed amongst non-zero land type fraction fields.

Fields returned:
- JULES land cover type fraction ancillary fields ('m01s00i216').
- Land sea mask ('m01s00i030')
- Land area fractions ('m01s00i505')

"""
import os
import warnings

import ants
import ants.decomposition as decomp
import ants.fileformats.cover_mapping as cover_mapping
import ants.io.save as save
import ants.utils
import iris
import numpy as np
from proc_ants import lct


def _load_ostia(ostia):
    # DEV: This function is here as part of development but isn't currently used.
    #      The OSTIA mask cannot be loaded without processing due to bad metadata.
    # Note that there may be at some point some cause to interpret the OSTIA
    # mask in some way and make some conditional on the lakes included in the
    # mask based on their presence in OSTIA.
    #
    # For example, we can fish out lat-lon points for each lake geometry in
    # Natural Earth as follows:
    # >>> records = proc_ants.lakes.get_lake_geoms()
    # >>> for record in records:
    #    geom = record.geometry.buffer(0)
    #    x, y = geom.representative_point().xy
    #    x = x[0]
    #    y = y[0]
    #
    # However, we should note that not all lakes have identified geometries in
    # Natural Earth or necessarily with one that is aligned to data derived
    # from the CCI.  This means that relating lakes between OSTIA and the CCI
    # might be rather difficult to achieve programatically.
    #
    fields = iris.fileformats.pp.load(ostia)

    def mutate(field):
        field.lbproc = 0
        return field

    tweaked_fields = (mutate(field) for field in fields)
    ostia_cubes = ants.fileformats.pp.pp2cubes(tweaked_fields)
    assert len(ostia_cubes) == 1, "More than one OSTIA cube?"
    ostia_cube = ostia_cubes[0]
    ostia_lsm = ostia_cube.copy(np.ma.getmaskarray(ostia_cube.data))
    return ostia_lsm


def _prepare_mask_cube(cube):
    """
    Sets the attributes needed for mask files, and converts datatype to integer.

    Operates in place.

    Attributes 'valid_min' and 'valid_max' are set to ensure netCDF result is
    correct for boolean data.

    STASH code is set to m01s00i030.

    Data type is set to integer.

    Parameters
    ----------
    cube : :class:`iris.cube.Cube`
        Mask cube to prepare.

    Returns
    -------
    : None
    Operates in place.

    """
    cube.attributes["valid_min"] = 0
    cube.attributes["valid_max"] = 1
    cube.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i030")
    cube.data = cube.data.astype("int8")


def load_data(
    source_path,
    transform_path,
    target_grid=None,
    target_landseamask=None,
    land_fraction_threshold=None,
):
    source = ants.load_cube(source_path, "land_cover_lccs")
    if target_grid:
        target_cube = ants.load_grid(target_grid)
    else:
        target_cube = ants.fileformats.load_landsea_mask(
            target_landseamask, land_fraction_threshold
        )

    target_cube = target_cube.extract(ants.ExtractConstraint(source, fix_period=True))
    try:
        source = source.extract(ants.ExtractConstraint(target_cube))
    except:
        target_cube = target_cube.intersection(longitude=(-180, 180))
        source = source.extract(ants.ExtractConstraint(target_cube))
    trans = cover_mapping.load_cover_mapper(transform_path)
    return source, target_cube, trans


def _get_transposed_source(source):
    # Return a transposed array view so that we can utilise broadcasting rules
    # with masks.
    pseudo_level = source.coord("pseudo_level")
    pdim = source.coord_dims(pseudo_level)[0]
    transpose_indx = list(range(source.ndim))
    transpose_indx.pop(pdim)
    transpose_indx += [pdim]

    return ants.utils.ndarray.transposed_view(source.data, transpose_indx)


def set_whole_fraction_ice(source, threshold, ice_tile_id):
    """
    Set ice fractions to [0, 1].

    Parameters
    ----------
    source : `iris.cube.Cube`
        Source land cover type fraction, with pseudo-level coordinate
        representing the classes.

    See Also
    --------
    :func:`~ants.analysis.cover_mapping.normalise_fractions`

    """
    # Ensure that all fractions add up to 1 before applying the threshold.
    # This function currently assumed this to be true as the crosswalk
    # calculation already normalises before returning a result.
    # normalise_fractions(source)

    pseudo_level = source.coord("pseudo_level")
    if pseudo_level.ndim != 1:
        msg = 'Expecting 1D "peudo_level" coordinate, got {}D'
        raise ValueError(msg.format(pseudo_level.ndim))

    slices = ants.analysis.cover_mapping.fetch_lct_slices(source, ice_tile_id)
    ice_data = source.data[slices]

    data = np.asarray(_get_transposed_source(source))

    # Set data above threshold.
    mask = ice_data >= threshold
    data[mask] = 0
    ice_data[mask] = 1

    # Set data below threshold - this isn't the inverse of the mask above as
    # we may have 'nan' values present which we don't want to lose.
    mask = ice_data < threshold
    ice_data[mask] = 0
    ants.analysis.cover_mapping.normalise_fractions(source)


def remove_non_glacial_ice(lct, ice_id, soil_id):
    """
    Replace isolated ice and replace with soil.

    - Identify ice seed points by those locations with 5x5 surrounding ice
      cells.
    - Find contiguous ice regions using these seed points.
    - Replace ice with soil for those not matching these contiguous regions.

    """
    ice_level = ants.analysis.cover_mapping.fetch_lct_slices(lct, ice_id)
    soil_level = ants.analysis.cover_mapping.fetch_lct_slices(lct, soil_id)
    ice_level_cube = lct[ice_level]

    # Find ice seed points (5x5 areas of ice)
    mask = np.ones((5, 5))
    ice_seed = ants.analysis.MooreNeighbourhood(
        ice_level_cube, window_mask=mask
    ).all_equal_value(1)

    # Use seed points (above) to find contiguous ice regions and construct a
    # mask which excludes these contiguous ice regions.
    ice_data = lct.data[ice_level]
    replace_ice_mask = ice_data > 0  # all ice assigned for removal.
    while ice_seed.any():
        y, x = np.unravel_index(ice_seed.argmax(), ice_seed.shape)
        ice_contiguous = ants.analysis.find_similar_region(
            ice_data,
            (y, x),
            extended_neighbourhood=True,
            wraparound=lct.coord(axis="x").circular,
        )
        # Remove contiguous region from seed array.
        ice_seed[ice_contiguous] = False
        # Remove contiguous ice from ice removal mask.
        replace_ice_mask[ice_contiguous] = False

    # Replace ice with soil.
    soil_data = lct.data[soil_level]
    soil_data[replace_ice_mask] = 1
    ice_data[replace_ice_mask] = 0


def gen_lct(src_cube, grid_cube, transform, min_frac, ice_tile_id, soil_tile_id):
    operation = ants.analysis.SCTTransformer(transform)
    lct_cube = decomp.decompose(operation, src_cube, grid_cube)
    lct_cube.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i216")
    lct_cube.rename("vegetation_area_fraction")

    # Any masked data is nan in the land cover type fraction as the ocean level
    # is used to derive the landsea mask.
    if np.ma.is_masked(lct_cube.data):
        lct_cube.data[lct_cube.data.mask] = np.nan
        lct_cube.data = lct_cube.data.data

    # Derive the land area fraction before we remove the ocean.
    # This field is not currently used directly by the model.
    land_fraction_cube = lct_cube.extract(iris.Constraint(pseudo_level=0))
    land_fraction_cube.remove_coord("pseudo_level")
    land_fraction_cube.data = 1 - land_fraction_cube.data
    land_fraction_cube.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi(
        "m01s00i505"
    )
    land_fraction_cube.rename("land_area_fraction")

    set_whole_fraction_ice(lct_cube, 0.5, ice_tile_id)
    remove_non_glacial_ice(lct_cube, ice_tile_id, soil_tile_id)
    # Exclude the ocean level - All land fractions must add up to 1 which means
    # that points must be 100% land or 100% ocean.  The array is masked as a
    # result (representing this ocean).
    lct_cube, land_mask_cube = lct.remove_ocean_level(lct_cube, min_frac=min_frac)

    # Extract landsea mask.
    land_mask_cube = lct_cube.slices_over("pseudo_level").next()
    land_mask_cube.remove_coord("pseudo_level")
    # Ensure land_mask_cube is appropriately masked
    ants.utils.cube.fix_mask(land_mask_cube)
    land_mask_cube.data = np.logical_not(land_mask_cube.data.mask).astype(
        "int8", copy=False
    )
    nan_values = np.isnan(lct_cube.data[0])
    # Inherit unknown values from the lct in the form of a mask.
    if nan_values.any():
        land_mask_cube.data = np.ma.array(land_mask_cube.data, mask=nan_values)
        warnings.warn(
            "There are nan values in the vegetation fraction. These will translate "
            "to masked data in the land sea mask and may cause unexpected results "
            "if the mask is used with ants.analysis.make_consistent_with_lsm."
        )

    _prepare_mask_cube(land_mask_cube)
    land_mask_cube.rename("land_binary_mask")
    lsm_cubes = iris.cube.CubeList([land_mask_cube, land_fraction_cube])
    return lct_cube, lsm_cubes

def main(
    source_path,
    transform_path,
    out_filepath,
    target_grid,
    landseamask_in,
    landseamask_out_root,
    land_fraction_threshold,
    ice_tile_id,
    soil_tile_id,
):
    source, grid, src_trans, = load_data(
        source_path,
        transform_path,
        target_grid,
        landseamask_in,
        land_fraction_threshold,
    )

    # Select land fraction to mask out
    min_frac = 0.1
    if landseamask_in:
        min_frac = 0.0

    lct_cube, lsm_cubes = gen_lct(source, grid, src_trans, min_frac, ice_tile_id, soil_tile_id)

    if landseamask_out_root:
        ants.config.dirpath_writeable(landseamask_out_root)
        _prepare_mask_cube(lsm_cubes[0])
        land_mask = lsm_cubes[0]
        land_fraction = lsm_cubes[1]
        ocean_data = np.ones_like(land_fraction.data)
        ocean_data[land_fraction.data < 1] = 0
        ocean_mask = land_fraction.copy(ocean_data)
        ocean_mask.standard_name = "land_binary_mask"
        _prepare_mask_cube(ocean_mask)

        cubes_to_save = {
            # "filename": (cube, fill_value),
            "qrparm.mask": (land_mask, -1),
            "qrparm.mask_sea": (ocean_mask, -1),
            "qrparm.landfrac": (land_fraction, None),
        }
        for filename, (cube, fill_value) in cubes_to_save.items():
            filepath = os.path.join(landseamask_out_root, filename)
            save_ants(cube, filepath, fill_value)

    if landseamask_in:
        ants.analysis.make_consistent_with_lsm(lct_cube, grid, True)
    save_ants(lct_cube, out_filepath, None)
    return lct_cube, lsm_cubes[0], lsm_cubes[1]


def save_ants(cube, filename, fill_value):
    save.ancil(cube, filename)
    save.netcdf(cube, filename, fill_value=fill_value)


def _get_parser():
    parser = ants.AntsArgParser(target_lsm=True, target_grid=True)
    parser.add_argument(
        "--transform-path", type=str, help="Filepath to crosswalk table.", required=True
    )
    msg = (
        "Output directory for the landseamask as derived by the the land cover type "
        "fraction field."
    )
    parser.add_argument("--landseamask-output-root", type=str, help=msg, required=False)
    parser.add_argument("--ice-tile-id", type=int, help="ID of the ice tile", default=9)
    parser.add_argument("--soil-tile-id", type=int, help="ID of the bare soil tile", default=8)
    return parser


if __name__ == "__main__":
    args = _get_parser().parse_args()
    if args.landseamask_output_root and args.target_lsm:
        raise ValueError(
            "Output landseamask filepath (landseamask_filepath) "
            "provided, along with an input landsea mask "
            "(target-lsm).  These are mutually exclusive "
            "arguments."
        )
    main(
        args.sources,
        args.transform_path,
        args.output,
        args.target_grid,
        args.target_lsm,
        args.landseamask_output_root,
        args.land_threshold,
        args.ice_tile_id,
        args.soil_tile_id,
    )
