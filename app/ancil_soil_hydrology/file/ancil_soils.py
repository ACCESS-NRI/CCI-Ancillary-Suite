#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
"""
Soils hydrology, carbon and thermal property parameter generation
*****************************************************************

- Calculate soils hydrology parameters by providing a soils field containing
  a geo-located field describing an 'MU_GLOBAL' at a given location.  We also
  provide a json soil lookup table which maps these MU_GLOBAL to unique soil
  IDs which are in turn mapped to various parameters (inc. sand-silt-clay
  fraction).  This lookup contains two table lookups:
- Define a unique soil ID for each unique soil sand-silt-clay fraction in the
  lookup.
    - muglobal_soilid_lookup: A dictionary mapping MU_GLOBAL to unique soil
      IDs along with the SHARE (percentage covergae) corresponding to this soil
      ID.
    - soilid_parameter_lookup: A dictionary mapping unique soil ID to the
      actual parameters in the database (sand, silt, clay etc. parameters as
      well as Cosby parameters).  This lookup and the exact values contained
      within it depend on the hydraulics scheme used to derive these Cosby
      parameters.

  To illustrate::

    {'muglobal_soilid_lookup': {mu_global_1: {unique_soil_id_1: share_1,
                                              ...,
                                              unique_soil_id_2: share_2},
                                mu_global_2: {unique_soil_id_1: share_1,
                                              ...,
                                              unique_soil_id_2: share_2}},
     'soilid_parameter_lookup': {unique_soil_id_1: {soil input parameters inc
                                                    sand, silt & clay},
                                 ...,
                                 unique_soil_id_n: {soil input parameters inc
                                                    sand, silt & clay}},}

- Determine the dominant sand-silt-clay fractions for each target cell.
- Determine the hydrology parameter values for each target cell using the
  lookup for Cosby parameters.
- Ensure consistency with the land cover type fraction (ice) ancillary:
    - Coastline adjustment.
    - Ice value override.
        - Set soil parameter values at locations of ice to 0 for all fields
          except `soil_thermal_capacity` and `soil_thermal_conductivity` which
          are set to 630000 and 0.2650 respectively.
- Save the soils hydrology parameter ('m01s00i207', 'm01s00i048',
  'm01s00i044', 'm01s00i043', 'm01s00i040', 'm01s00i041', 'm01s00i008'),
  thermal properties ('m01s00i047', 'm01s00i046') and soil carbon
  ('m01s00i223') cubes.
- Additionaly saves volume fraction of sand/silt/clay in soil ('m01s00i420',
  'm01s00i419' and 'm01s00i418' respectively).  Not directly used by the model.

"""
import functools
import json
from abc import ABCMeta, abstractmethod

import iris
import numpy as np

import ants
import ants.decomposition as decomp
import ants.io.save as save
import ants.utils
import soil_cosby_parameters

_ICE_VALUE = {"soil_thermal_capacity": 630000, "soil_thermal_conductivity": 0.2650}


# jules parameters to be saved to the ancillary
MAPPING = {
    "b": {
        "units": 1,
        "STASH": "m01s00i207",
        "standard_name": None,
        "long_name": ("Brooks-Corey exponent in the soil hydraulic " "characteristics"),
    },
    "sathh": {
        "units": "m",
        "STASH": "m01s00i048",
        "standard_name": "soil_suction_at_saturation",
        "long_name": ("absolute value of the soil matric suction at " "saturation"),
    },
    "satcon": {
        "units": "kg m-2 s-1",
        "STASH": "m01s00i044",
        "standard_name": "soil_hydraulic_conductivity_at_saturation",
        "long_name": "hydraulic conductivity at saturation",
    },
    "sm_sat": {
        "units": "m3 m-3",
        "STASH": "m01s00i043",
        "standard_name": "soil_porosity",
        "long_name": "volumetric soil moisture content at saturation",
    },
    "sm_wilt": {
        "units": "m3 m-3",
        "STASH": "m01s00i040",
        "standard_name": (
            "volume_fraction_of_condensed_water_in_soil_" "at_wilting_point"
        ),
        "long_name": ("volumetric soil moisture content at the " "wilting point"),
    },
    "sm_crit": {
        "units": "m3 m-3",
        "STASH": "m01s00i041",
        "standard_name": (
            "volume_fraction_of_condensed_water_in_soil_" "at_critical_point"
        ),
        "long_name": (
            "volumetric soil moisture content at the "
            "critical point at which soil moisture stress "
            "starts to reduce transpiration"
        ),
    },
    "BD": {
        "units": "kg m-3",
        "STASH": "m01s00i008",
        "long_name": ("bulk density of the soil"),
    },
    "soil_carb": {
        "units": "kg m-2",
        "STASH": "m01s00i223",
        "standard_name": "soil_carbon_content",
        "long_name": ("amount of soil carbon in the soil"),
    },
    "hcon": {
        "units": "W m-1 K-1",
        "STASH": "m01s00i047",
        "standard_name": "soil_thermal_conductivity",
        "long_name": "dry thermal conductivity (lambda)",
    },
    "hcap": {
        "units": "J m-3 K-1",
        "STASH": "m01s00i046",
        "standard_name": "soil_thermal_capacity",
        "long_name": "dry heat capacity",
    },
}


# Not strictly jules parameters but to be saved to the ancillary
MAPPING_EXTENDED = MAPPING.copy()
MAPPING_EXTENDED.update({
    'T_SAND': {'standard_name': 'volume_fraction_of_sand_in_soil',
               'long_name': 'volume_fraction_of_sand_in_soil',
               'units': 'm3 m-3', 'STASH': 'm01s00i420'},
    'T_SILT': {'standard_name': 'volume_fraction_of_silt_in_soil',
               'long_name': 'volume_fraction_of_silt_in_soil',
               'units': 'm3 m-3', 'STASH': 'm01s00i419'},
    'T_CLAY': {'standard_name': 'volume_fraction_of_clay_in_soil',
               'long_name': 'volume_fraction_of_clay_in_soil',
               'units': 'm3 m-3', 'STASH': 'm01s00i418'}
})

class _SoilRegridScheme(object, metaclass=ABCMeta):
    """Abstract class which defined the UI of the Soil regrid scheme."""

    @abstractmethod
    def _ret_parameters(self, lookup, units_cube, target_grid):
        """
        Calculate soils hydrology, carbon and thermal property parameters.

        Parameters
        ----------
        lookup : :class:`numpy.ndarray`
            Soil lookup table.
        units_cube : :class:`~iris.cube.Cube`
            Soils MU_GLOBAL field.
        target_grid : :class:`~iris.cube.Cube`
            Definition of the output target grid.

        Returns
        -------
        : iterable of :class:`~iris.cube.Cube` objects.
            Soils hydrology, carbon and thermal property parameters.

        """
        raise NotDefinedError

    @classmethod
    def __call__(cls, lookup, units_cube, target_grid):
        """
        Calculate soils hydrology, carbon and thermal property parameters.

        Parameters
        ----------
        lookup : :class:`numpy.ndarray`
            Soil lookup table.
        units_cube : :class:`~iris.cube.Cube`
            Soils MU_GLOBAL field.
        target_grid : :class:`~iris.cube.Cube`
            Definition of the output target grid.

        Returns
        -------
        : iterable of :class:`~iris.cube.Cube` objects.
            Soils hydrology, carbon and thermal property parameters.

        """
        return cls._ret_parameters(lookup, units_cube, target_grid)

    @staticmethod
    def _weights(source, target):
        """
        Return ESMF weights for the given source-target pair.

        Parameters
        ----------
        source : :class:`~iris.cube.Cube`
            Source dataset to regrid from.
        target : :class:`~iris.cube.Cube`
            Target dataset to regrid to.

        Returns
        -------
        : `numpy.ndarray`, `numpy.ndarray`, `numpy.ndarray`
            source indices (columns), target indices (rows) and regridding
            weights between source and target cells.

        """
        scheme = ants.regrid.esmf.ConservativeESMF()
        regridder = scheme.regridder(source, target, persistent_cache=True)
        columns, rows, weights = regridder.cache
        return columns, rows, weights


class DomParameters(_SoilRegridScheme):
    """
    Return dominant soil parameters for each target cell.

    Calculate the area weights between the MU_GLOBAL field and the target grid,
    then lookup which soil IDs contribute to a given target grid cell.
    From this, the share for each of these soil ID's is multiplied by the area
    weight.  The soil ID with the greatest area weight x share is then deemed
    the dominant soil type.  The parameters returned then result from a lookup
    in the table using this soil ID.

    """

    @classmethod
    def _fetch_param_cubes(cls, lookup, dom_soil_id):
        """
        Return the soils hydrology, carbon and thermal property parameter set
        from the dominant soil ID field.

        For a soil ID at a given target cell, use the lookup to return all
        relevant soils hydrology, carbon and thermal property parameters as a
        set of cubes.

        Parameters
        ----------
        lookup : :class:`numpy.ndarray`
            Soil lookup table.
        dom_soil_id : :class:`~iris.cube.Cube`
            A field containing a unique soil ID at each location.

        Returns
        -------
        : iterable of :class:`~iris.cube.Cube` objects.
            Soils hydrology, carbon and thermal property parameters.

        """
        soilid_parameter_lookup = lookup["soilid_parameter_lookup"]

        cubes = iris.cube.CubeList(
            [
                dom_soil_id.copy(np.ma.zeros(dom_soil_id.shape))
                for i in range(len(list(MAPPING_EXTENDED.keys())))
            ]
        )
        cubes_dict = {}
        for name, cube in zip(list(MAPPING_EXTENDED.keys()), cubes):
            try:
                cube.standard_name = MAPPING_EXTENDED[name]['standard_name']
            except KeyError:
                pass
            cube.long_name = MAPPING_EXTENDED[name]['long_name']
            cube.var_name = name
            cube.units = MAPPING_EXTENDED[name]['units']
            cube.attributes['STASH'] = iris.fileformats.pp.STASH.from_msi(
                MAPPING_EXTENDED[name]['STASH'])
            cubes_dict[name] = cube

        for i in range(dom_soil_id.shape[0]):
            for j in range(dom_soil_id.shape[1]):
                soil_id = dom_soil_id.data[i, j]
                if soil_id in [0, np.ma.masked]:
                    continue
                params = soilid_parameter_lookup[str(soil_id)]
                for key, value in params.items():
                    try:
                        cubes_dict[key].data[i, j] = value
                    except KeyError:
                        pass
        if np.ma.isMaskedArray(dom_soil_id.data):
            for cube in cubes:
                cube.data.mask = dom_soil_id.data.mask.copy()
        return cubes

    @staticmethod
    def _ret_dominant_units(lookup, units, rows, columns, weights, shape):
        """
        Determine the dominant soil ID for each target cell.

        Parameters
        ----------
        lookup : :class:`numpy.ndarray`
            Soil lookup table.
        units : :class:`~iris.cube.Cube`
            Soils MU_GLOBAL field.
        rows : `numpy.ndarray`
            source indices (columns).
        columns : `numpy.ndarray`
            target indices (rows)
        weights : `numpy.ndarray`
            Regridding weights between source and target cells.
        shape : tuple
            Target grid shape.

        Returns
        -------
        : :class:`numpy.ndarray`
            Dominant soil ID at each target grid cell location.

        """
        units_flattened = units.reshape(-1)
        muglobal_soilid_lookup = lookup["muglobal_soilid_lookup"]
        output = np.ma.zeros(shape, dtype="int")
        output.mask = True
        output_flattened = output.reshape(-1)

        curr_row = -1
        contributing = {}
        for row, column, weight in zip(rows, columns, weights):
            # Reset the maximum weight found for this target cell and move on
            # to recording the maximum weight for the next target cell.
            if row != curr_row:
                if contributing:
                    output_flattened[curr_row] = list(contributing.keys())[
                        np.array(list(contributing.values())).argmax()
                    ]
                    contributing.clear()
                curr_row = row

            # Fetch the MU_GLOBAL for this source cell.
            mu_global = units_flattened[column]
            if not mu_global > 0:
                # skip ocean
                continue
            # Fetch all soil IDs corresponding to this MU_GLOBAL
            munits = muglobal_soilid_lookup[str(mu_global)]
            for unit in munits:
                if int(unit) == 1:
                    # skip non-soil
                    continue
                # Determine contribution of this soil ID (weight * share).
                if unit not in contributing:
                    contributing[unit] = 0
                contributing[unit] += weight * munits[unit]

        if contributing:
            output_flattened[curr_row] = list(contributing.keys())[
                np.array(list(contributing.values())).argmax()
            ]
        return output

    @classmethod
    def _ret_parameters(cls, lookup, units_cube, target_grid):
        """
        Calculate soils hydrology parameters.

        Parameters
        ----------
        lookup : :class:`numpy.ndarray`
            Soil lookup table.
        units_cube : :class:`~iris.cube.Cube`
            Soils MU_GLOBAL field.
        target_grid : :class:`~iris.cube.Cube`
            Definition of the output target grid.

        Returns
        -------
        : iterable of :class:`~iris.cube.Cube` objects.
            Soils hydrology, carbon and thermal property parameters.

        """
        columns, rows, weights = cls._weights(units_cube, target_grid)
        dom_units = cls._ret_dominant_units(
            lookup, units_cube.data, rows, columns, weights, target_grid.shape
        )
        dom_units_cube = target_grid.copy(dom_units)
        return cls._fetch_param_cubes(lookup, dom_units_cube)


class AvgParameters(_SoilRegridScheme):
    """
    Return the average soil parameters for each target cell.

    Calculate the area weighted mean sand, silt, clay, bulk density and
    organic chemistry values for a given target cell.  Calculate the soil
    parameters using the relevant hydraulic scheme by passing these values as
    input.

    """

    SOIL_HYDRAULICS_SCHEME = soil_cosby_parameters.ClappHornberger

    @classmethod
    def _fetch_param_cubes(cls, avg_params, target_cube):
        """
        Calculate cosby parameters for the provided average sand, silt, clay,
        bulk density and organic carbon fraction parameters.

        Parameters
        ----------
        avg_params : :class:`numpy.ndarray`
            Average parameter set of shape (5, target_y, target_x), where '5'
            comes from the set of input parameters sand, silt, clay, bulk
            density and organic carbon fraction.
        target_cube : :class:`~iris.cube.Cube`
            Target grid definition.

        Returns
        -------
        : iterable of :class:`~iris.cube.Cube` objects.
            Soils hydrology, carbon and thermal property parameters.

        """
        cubes = iris.cube.CubeList(
            [
                target_cube.copy(np.ma.zeros(target_cube.shape))
                for i in range(len(list(MAPPING_EXTENDED.keys())))
            ]
        )
        cubes_dict = {}
        for name, cube in zip(list(MAPPING_EXTENDED.keys()), cubes):
            try:
                cube.standard_name = MAPPING_EXTENDED[name]['standard_name']
            except KeyError:
                pass
            cube.long_name = MAPPING_EXTENDED[name]['long_name']
            cube.var_name = name
            cube.units = MAPPING_EXTENDED[name]['units']
            cube.attributes['STASH'] = iris.fileformats.pp.STASH.from_msi(
                MAPPING_EXTENDED[name]['STASH'])
            cubes_dict[name] = cube

        for j in range(avg_params.shape[1]):
            for i in range(avg_params.shape[2]):
                if sum(avg_params[:3, j, i]) <= 0 or np.ma.is_masked(
                    avg_params[0, j, i]
                ):
                    for cube in cubes:
                        cube.data[j, i] = np.ma.masked
                    continue

                params = cls.SOIL_HYDRAULICS_SCHEME(
                    f_sand=avg_params[0, j, i],
                    f_silt=avg_params[1, j, i],
                    f_clay=avg_params[2, j, i],
                    bulk_den=avg_params[3, j, i],
                    org_carb=avg_params[4, j, i],
                )
                for name in list(MAPPING.keys()):
                    cubes_dict[name].data[j, i] = params[name]
        return cubes

    @classmethod
    def _ret_parameters(cls, lookup, units_cube, target_grid):
        """
        Calculate soils hydrology parameters.

        Parameters
        ----------
        lookup : :class:`numpy.ndarray`
            Soil lookup table.
        units_cube : :class:`~iris.cube.Cube`
            Soils MU_GLOBAL field.
        target_grid : :class:`~iris.cube.Cube`
            Definition of the output target grid.

        Returns
        -------
        : iterable of :class:`~iris.cube.Cube` objects.
            Soils hydrology, carbon and thermal property parameters.

        """
        columns, rows, weights = cls._weights(units_cube, target_grid)
        ret_params = ["T_SAND", "T_SILT", "T_CLAY", "T_BULK_DENSITY", "T_OC"]
        avg_params = cls._ret_average_params(
            lookup,
            units_cube.data,
            rows,
            columns,
            weights,
            target_grid.shape,
            ret_params,
        )
        return cls._fetch_param_cubes(avg_params, target_grid)

    @staticmethod
    def _ret_average_params(lookup, units, rows, columns, weights, shape, ret_params):
        """
        Determine the average parameters for each target cell.

        Iterate over each source cell contributing over each target cell and
        calculate the average fractions requested in 'ret_params'.

        Parameters
        ----------
        lookup : :class:`numpy.ndarray`
            Soil lookup table.
        units : :class:`~iris.cube.Cube`
            Soils MU_GLOBAL field.
        rows : `numpy.ndarray`
            source indices (columns).
        columns : `numpy.ndarray`
            target indices (rows)
        weights : `numpy.ndarray`
            Regridding weights between source and target cells.
        shape : tuple
            Target grid shape.
        ret_params : list
            List of parameter names used to derive the soils parameter set.

        Returns
        -------
        : :class:`numpy.ndarray`
            Average parameter set of shape (5, target_y, target_x), where '5'
            comes from the set of input parameters sand, silt, clay, bulk
            density and organic carbon fraction.

        """

        def _ensure_add_up_to_1(sand_silt_clay):
            """Ensure the T_SAND, T_SILT and T_CLAY fractions add to 1."""
            # Adjust to ensure they add up to 1 (add difference to dominant
            # type).
            diff = 1.0 - sum(sand_silt_clay)
            if diff > 0.1:
                msg = "sand-silt-clay sum does not add up to 1: {}"
                raise ValueError(msg.format(diff))
            argmax = sand_silt_clay.argmax()
            sand_silt_clay[argmax] += diff

        units_flattened = units.reshape(-1)
        muglobal_soilid_lookup = lookup["muglobal_soilid_lookup"]
        soilid_parameter_lookup = lookup["soilid_parameter_lookup"]

        # Keep only those average params. used to derive the output cosby set.
        headers = list(list(soilid_parameter_lookup.values())[0].keys())
        headers = [headers.index(ret_param) for ret_param in ret_params]

        output = np.ma.zeros([len(headers)] + list(shape), dtype="float")
        output.mask = True
        # We flatten everything except for the parameter dimension.
        output_flattened = output.reshape(output.shape[0], -1)

        curr_row = -1
        total_weight = 0
        for row, column, weight in zip(rows, columns, weights):

            # Reset the total weight found if we have for this particular
            # target cell (curr_row).
            # Calculate the average parameters.
            # Continue to collect weights for next target cell (curr_row).
            if row != curr_row:
                if total_weight > 0:
                    output_flattened[:, curr_row] = (
                        output_flattened[:, curr_row] / total_weight
                    )
                    _ensure_add_up_to_1(output_flattened[:3, curr_row])
                    total_weight = 0
                curr_row = row

            # Fetch the MU_GLOBAL for this source cell.
            mu_global = units_flattened[column]
            if not int(mu_global) > 0:
                # skip ocean
                continue
            # Fetch all soil IDs corresponding to this MU_GLOBAL
            units = muglobal_soilid_lookup[str(mu_global)]
            for unit in units:
                if int(unit) == 1:
                    # skip non-soil
                    continue
                # Determine contribution of this soil ID (weight * share).
                sweight = weight * units[unit]
                total_weight += sweight
                params = np.array(list(soilid_parameter_lookup[unit].values()))
                output_flattened.mask[:, row] = False
                output_flattened[:, row] += sweight * params[headers]
        if total_weight > 0:
            output_flattened[:, curr_row] = output_flattened[:, curr_row] / total_weight
            _ensure_add_up_to_1(output_flattened[:3, curr_row])
        return output


def make_consistent_lct(parameters, lct_fractions, ice_tile_id=9):
    """
    Ensure consistency with the ice from the land cover types.

    - Mask locations where all parameters are zero at a given location.
    - Mask parameters fields where the lct_fractions identifies locations of
      ice.
    - Use the spiral search in making the mask consistent with the mask from
      the provided lct_fractions.
    - Override ice values to 0 in all parameter fields except
      `soil_thermal_capacity` and `soil_thermal_conductivity` which are set to
      630000 and 0.2650 respectively.

    Parameters
    ----------
    parameters : iterable of :class:`~iris.cube.Cube` objects.
        Each of the Cosby parameter set.
    lct_fractions : :class:`~iris.cube.Cube`
        Land cover fraction fields, which amongst these include the fraction
        of ice cover as well as the mask on which this applies.

    Returns
    -------
    : None
        In-place operation

    """
    # Override parameter ice values with 0 to ensure Fill doesn't fill them.
    ice_slices = ants.analysis.cover_mapping.fetch_lct_slices(
        lct_fractions, ice_tile_id
    )
    ice_level_cube = lct_fractions[ice_slices]

    # Mask locations where the all parameters are 0 at that location.
    all_param_non_zero = np.zeros(parameters[0].shape, dtype="bool")
    for param in parameters:
        all_param_non_zero += param.data != 0
    all_zero = ~all_param_non_zero
    for param in parameters:
        param.data[all_zero] = np.ma.masked

    # Ensure ice_level_cube is masked
    ants.utils.cube.fix_mask(ice_level_cube)

    # lsm + ice
    lsm = ice_level_cube.copy(np.ma.getmaskarray(ice_level_cube.data))
    ice_mask = (ice_level_cube.data.data > 0) & (~ice_level_cube.data.mask)
    target_mask_noice = lsm.copy(ice_mask + lsm.data)

    # Make consistent with land cover type fraction landsea mask.
    filler = ants.analysis.FillMissingPoints(
        parameters[0], target_mask=target_mask_noice
    )
    for param in parameters:
        filler(param)
        if param.name() in _ICE_VALUE:
            # Override ICE value.
            param.data[ice_mask] = _ICE_VALUE[param.name()]
        else:
            param.data[ice_mask] = 0


def load(units_filepath, lct_fraction_filepath):
    """
    Load the input datasets.

    Parameters
    ----------
    units_filepath : str
        Input soils dataset defining the geo-located MU_GLOBAL at each
        location.
    lct_fraction_filepath : str
        Vegetation fraction ancillary filepath.

    Results
    -------
    : :class:`~iris.cube.Cubes`, :class:`~iris.cube.Cubes`
        Return both soils MU_GLOBAL and land covert type fraction dataset.

    """
    units_cube = ants.io.load.load_cube(units_filepath)

    stash_con = iris.AttributeConstraint(STASH="m01s00i216")
    lct_fractions = ants.io.load.load_cube(lct_fraction_filepath, stash_con)
    ants.utils.cube.guess_horizontal_bounds(lct_fractions)
    return units_cube, lct_fractions


def soils_gen(lookup_table, units_cube, target_grid):
    """
    Generate soils parameters.

    Load the soils lookup table and determine the dominate soil type at each
    target grid cell.

    Parameters
    ----------
    lookup_table : :class:`numpy.ndarray`
        Soils lookup table.
    units_cube : :class:`~iris.cube.Cube`
        Soils MU_GLOBAL field.
    target_grid : :class:`~iris.cube.Cube`
        Definition of the output target grid.

    Returns
    -------
    : iterable of :class:`~iris.cube.Cube` objects.
        Soils hydrology, carbon and thermal property parameters.

    """
    lookup_table = lookup_table_loader(lookup_table)
    return DomParameters()(lookup_table, units_cube, target_grid)


def lookup_table_loader(lookup_filepath):
    """
    Lookup table loader.

    A loader is passed to the operation to minimise data transfer between the
    operating process and the controlling process.

    Parameters
    ----------
    lookup_filepath : str
        Filepath to the input soils lookup table.

    Return
    ------
    : Column named `numpy.ndarray`
        Soil lookup table.

    """
    with open(lookup_filepath, "r") as fh:
        lookup_table = json.load(fh)
    return lookup_table


def main(soil_units_filepath, lct_fraction_filepath, lookup_filepath, ice_tile_id, output, netcdf_only):
    """
    Main calling function to calculate the soil parameters.

    Calculate the dominant Cosby parameters for each target cell location and
    ensure consistency with the pre-existing ice and landsea mask fields.

    Parameters
    ----------
    soils_units_filepath : str
        Input soils dataset filepath defining the geo-located MU_GLOBAL at each
        location.
    lct_fraction_filepath : str
        Vegetation fraction ancillary filepath.
    lookup_filepath : str
        Filepath to the input soils lookup table.

    Returns
    -------
    : None

    """
    units_cube, lct_fractions = load(soil_units_filepath, lct_fraction_filepath)

    operation = functools.partial(soils_gen, lookup_filepath)
    grid = lct_fractions[0]
    grid.remove_coord("pseudo_level")
    parameters = decomp.decompose(operation, units_cube, grid)

    # Make the unique soil ID cube consistent with the land cover type
    # fractions (lct_fractions) ancillary.
    make_consistent_lct(parameters, lct_fractions, ice_tile_id)

    save.netcdf(parameters, output)
    if not netcdf_only:
        save.ancil(parameters, output)


if __name__ == "__main__":
    parser = ants.AntsArgParser()
    parser.add_argument(
        "--lct-ancillary",
        type=str,
        help="Filepath to the land covert type fraction " "ancillary.",
        required=True,
    )
    parser.add_argument(
        "--soils-lookup",
        type=str,
        help="Soils hydrology lookup between source mapping "
        " and corresponding hwsd parameters.",
        required=True,
    )
    parser.add_argument(
        "--ice-tile-id",
        type=int,
        help="ID for the ice tile",
        default=9
    )
    args = parser.parse_args()
    main(args.sources, args.lct_ancillary, args.soils_lookup, args.ice_tile_id, args.output, args.netcdf_only)
