# This creates basic soil initial conditions, that should be valid enough
# to begin a spin-up of a land model. Should not be used for any short time
# scale simulations, where the initial condition is important.

import argparse
import numpy

import iris


def _parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
            '--soil-thickness',
            required=True,
            type=str,
            help='Comma separated list of soil layer thicknesses.'
            )

    parser.add_argument(
            '--soil-parameters',
            required=True,
            type=str,
            help='Path to the soil parameters file.'
            )

    parser.add_argument(
            '--moisture-units',
            required=False,
            type=str,
            default='volumetric',
            help="""Whether to output the soil moisture in 'volumetric' or
            'mass' units."""
            )

    parser.add_argument(
            '--max-temperature',
            required=False,
            type=float,
            default=25.0,
            help='Maximum soil temperature to apply at the equator.'
            )

    parser.add_argument(
            '--output',
            required=True,
            type=str,
            help='Path to write the output.'
            )

    return parser.parse_args()


def depths_from_thickness(soil_thickness):
    """
    Generate the soil layer depths from the given soil layer thickness.
    """
    return [0.5 * soil_thickness[l] if l == 0 else
              numpy.sum(soil_thickness[:l]) + 0.5 * soil_thickness[l]
              for l in range(len(soil_thickness))
              ]


def generate_soil_moisture(
        vol_smc_at_saturation,
        soil_thickness,
        moisture_units
        ):
    """
    Create a new soil moisture cube from the volumetric soil moisture
    content at saturation, given soil depths and desired moisture units. The
    resulting soil moisture is equivalent to the saturation fraction.

    Returns a new cube of "soil_moisture_content_in_a_layer".
    """

    # Broadcast the single dimensional soil properties over layers
    smc_at_sat = numpy.repeat(
            vol_smc_at_saturation.data[:, :, numpy.newaxis],
            len(soil_thickness),
            axis=2
            )

    if moisture_units == 'mass':
        # We want soil moisture in kg m-2, so scale by 
        # density * soil layer depth.
        water_density = 1000.0
        smc_at_sat = smc_at_sat * soil_thickness * water_density

        standard_name = 'mass_content_of_water_in_soil_layer'
        units = 'kg m-2'

    elif moisture_units == 'volumetric':
        standard_name = 'volumetric_soil_moisture_content_in_a_layer'
        units = 'm3 m-3'

    else:
        raise ValueError("""The supplied moisture units must be either
        'mass' or 'volumetric'.""")

    depths = depths_from_thickness(soil_thickness)

    # Create the soil depth dimension
    depth_coord = iris.coords.DimCoord(
            depths,
            standard_name='depth',
            units='m'
            )

    smc_cube = iris.cube.Cube(
            smc_at_sat,
            var_name=standard_name,
            units=units,
            dim_coords_and_dims=[
                (vol_smc_at_saturation.coord('latitude'), 0),
                (vol_smc_at_saturation.coord('longitude'), 1),
                (depth_coord, 2)
                ]
            )

    return smc_cube


def generate_soil_temperature(vol_smc_at_saturation, soil_depths, max_temp):
    """Generate a new soil temperature cube using a simple latitude dependence.
    The soil temperature varies from the supplied max temperature at the
    equator to 0.0 at the poles."""

    # Reuse the vol_smc_at_saturation so that we have the mask
    soil_temperature = numpy.ones_like(vol_smc_at_saturation.data)

    # Set up the latitude scaling
    latitude_scaling = (90 - numpy.abs(
            vol_smc_at_saturation.coord('latitude').points
    )) / 90

    # Apply the scaling and the max temperature
    soil_temperature = soil_temperature * \
            numpy.reshape(latitude_scaling, (-1, 1)) * max_temp

    # Spread over the soil layers
    soil_temperature = numpy.repeat(
            soil_temperature[:, :, numpy.newaxis],
            len(soil_thickness),
            axis=2
            )

    depths = depths_from_thickness(soil_thickness)

    # Create the soil depth dimension
    depth_coord = iris.coords.DimCoord(
            depths,
            standard_name='depth',
            units='m'
            )

    soil_temp_cube = iris.cube.Cube(
        soil_temperature,
        var_name='soil_temperature',
        units='K',
        dim_coords_and_dims=[
            (vol_smc_at_saturation.coord('latitude'), 0),
            (vol_smc_at_saturation.coord('longitude'), 1),
            (depth_coord, 2)
            ]
        )

    return soil_temp_cube


def main(soil_parameters_path,
         soil_thickness,
         moisture_units,
         max_temp,
         output
         ):
    # Load in the smc wilting to determine the mask and create smc
    smc_at_wilting = iris.load_cube(
            soil_parameters_path,
            'volumetric_soil_moisture_content_at_saturation'
            )

    soil_smc = generate_soil_moisture(
            smc_at_wilting,
            soil_thickness,
            moisture_units
            )

    soil_temp = generate_soil_temperature(
            smc_at_wilting,
            soil_thickness,
            max_temp
            )

    cubelist = iris.cube.CubeList([soil_smc, soil_temp])

    iris.fileformats.netcdf.save(cubelist, output)


if __name__ == '__main__':

    args = _parse_args()

    soil_thickness = numpy.array([float(d) 
                                  for d in args.soil_thickness.split(',')])

    main(
            args.soil_parameters,
            soil_thickness,
            args.moisture_units,
            args.max_temperature,
            args.output
            )
