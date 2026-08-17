import argparse
import ants
import iris
import numpy
import scipy

R_gas = 287.058
grav = 9.81
sea_level_pressure = 1013.0
mass_ozone = 48.0
mass_air = 28.97
ppmv2mmr = 1.65

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--ozone-path',
            required=True,
            type=str,
            help='Path to the ozone source. Can be a directory or file.'
            )
    parser.add_argument(
            '--temperature-path',
            required=True,
            type=str,
            help='Path to the temperature source. Can be a directory or file.'
            )
    parser.add_argument(
            '--orography-file',
            required=True,
            type=str,
            help='Path to the orography file containing surface_altitude.'
            )
    parser.add_argument(
            '--vertical-disc',
            required=True,
            type=str,
            help='Path to the file describing the vertical discretization.'
            )
    parser.add_argument(
            '--output',
            required=True,
            type=str,
            help='Path to write the output to.'
            )
    parser.add_argument(
            '--begin',
            required=False,
            type=int,
            help='Begin year for the output.'
            )
    parser.add_argument(
            '--end',
            required=False,
            type=int,
            help='End year for the output.'
            )
    parser.add_argument(
            '--zonal',
            action='store_true',
            help='Whether to reduce the output to zonal values.'
            )
    parser.add_argument(
            '--climatological',
            action='store_true',
            help='Whether to produce climatological output.'
            )
    parser.add_argument(
            '--extend',
            action='store_true',
            help="""Whether to allow constant extrapolation outside
            source data bounds."""
            )
    parser.add_argument(
            '--netcdf-only',
            action='store_true',
            help='Whether to only write to NetCDF, or include PP format as well.'
            )

    return parser.parse_args()


def main(
        ozone_path,
        temperature_path,
        orography_file,
        vert_file,
        output_file,
        begin_year,
        end_year,
        zonal,
        climatological,
        extend,
        netcdf_only
        ):

    ozone = iris.load(ozone_path)
    iris.util.equalise_attributes(ozone)
    ozone = ozone.concatenate_cube()

    temperature = iris.load(temperature_path)
    iris.util.equalise_attributes(temperature)
    temperature = temperature.concatenate_cube()

    orog_height = iris.load_cube(orography_file, 'surface_altitude')
    vert_nml = ants.fileformats.namelist._read_namelist(vert_file)
    vert_grid = ants.fileformats.namelist.VerticalLevels(vert_nml)

    new_ozone = generate_ozone(
            ozone,
            temperature,
            orog_height,
            vert_grid,
            begin_year,
            end_year,
            zonal,
            climatological,
            extend
            )

    ants.io.save.netcdf(new_ozone, output_file)
    if not netcdf_only:
        ants.io.save.ancil(new_ozone, output_file)

def generate_ozone(
        ozone,
        temperature,
        orog_height,
        vert_grid,
        begin_year,
        end_year,
        zonal=False,
        climatological=False,
        extend=False
        ):
    """
    Process the ozone data to produce a new dataset defined on model levels.
    """

    # Check that the supplied arguments are valid
    check_dates(ozone, begin_year, end_year)
    check_dates(temperature, begin_year, end_year)

    # Extract the desired years for ozone and optionally temperature
    yr_constraint = iris.Constraint(
                    time=lambda t: \
                    iris.time.PartialDateTime(begin_year) <= \
                    t.point <= \
                    iris.time.PartialDateTime(end_year)
                    )

    ozone = ozone.extract(yr_constraint)
    temperature = temperature.extract(yr_constraint)

    # Set up the real model level heights
    pressure_heights = derive_pressure_heights(temperature)
    pressure_heights_on_grid = interpolate_to_new_latitudes(
            pressure_heights,
            orog_height
            )
    model_heights = calculate_model_heights(orog_height, vert_grid)

    # Now interpolate the ozone onto the model horizontal grid
    match_coord_system(orog_height, ozone)
    ozone_on_horiz_grid = ozone.regrid(orog_height, iris.analysis.Linear())
    ozone_on_grid = vertical_interpolate(
                ozone_on_horiz_grid,
                pressure_heights_on_grid,
                model_heights
                )

    # Convert from molar fraction to mass fraction
    ozone_on_grid = ozone_on_grid * (mass_ozone / mass_air)

    # Perform zonal averaging if requested
    if zonal:
        ozone_on_grid = ozone_on_grid.collapsed(
                'longitude',
                iris.analysis.MEAN
                )

    # Perform climatological reduction. Note that this specifically averages
    # over the year range requested. If climatological is desired over the full
    # epoch of the original dataset, remove the --begin and --end arguments.
    if climatological:
        ozone_on_grid = iris.coord_categorisation.add_month_number(
                ozone_on_grid, 'time', name='month'
                )

        ozone_on_grid = ozone_on_grid.aggregation_by(
                'month',
                iris.analysis.MEAN,
                climatological=True
                )

    if extend:
        raise NotImplementedError("""Extending the ozone outside the range of
            the supplied data is not yet supported.""")

    return ozone_on_grid


def check_dates(cube, begin_year, end_year, extend=False):
    """
    Ensure that the data on the cube spans the given year range, assuming
    extend is False. If extend is True, just throw a message saying
    that the data is being extendd.
    """
    t_as_datetime = cube.coord(
            'time'
            ).units.num2date(cube.coord('time').points)
    yr_min = t_as_datetime[0].year
    yr_max = t_as_datetime[-1].year
    if yr_min > begin_year or yr_max < end_year:
        raise ValueError(f"""The cube's data spans from {yr_min} to 
        {yr_max}, but the range of years requested is {begin_year} to
        {end_year}.""")


def derive_pressure_heights(temperature):
    """
    Use the provided temperature cube, which has vertical coordinates of
    pressure levels, to turn them into heights.
    """
    
    pressure_levels = temperature.coord('pressure').points
    
    npres = len(pressure_levels)
    nt = len(temperature.coord('time').points)
    nlat = len(temperature.coord('latitude').points)

    integrant = R_gas * temperature / grav

    pressure_boundary = numpy.zeros(npres+1)
    pressure_boundary[0] = sea_level_pressure
    pressure_boundary[-1] = 1.0e-5
    pressure_boundary[1:-1] = numpy.sqrt(pressure_levels[:-1] * \
            pressure_levels[1:])

    log_pres_boundary = numpy.log(pressure_boundary)
    dlogpressure = numpy.log(pressure_boundary[:-1] / pressure_boundary[1:])
    frac = (log_pres_boundary[:-1] - numpy.log(pressure_levels)) / \
            (log_pres_boundary[:-1] - log_pres_boundary[1:])

    edge_height = numpy.zeros((nt, npres+1, nlat))
    centre_height = temperature.copy()

    for i in range(npres):
        edge_height[:, i+1, :] = edge_height[:, i, :] + dlogpressure[i] * \
                integrant.data[:, i, :]
        centre_height.data[:, i, :] = edge_height[:, i, :] + frac[i] * \
                dlogpressure[i] * integrant.data[:, i, :]

    return centre_height


def interpolate_to_new_latitudes(from_cube, to_cube):
    itp = iris.analysis.Linear(
            extrapolation_mode='linear'
            ).interpolator(from_cube, ['latitude'])
    interpolated = itp([to_cube.coord('latitude').points])

    return interpolated


def calculate_model_heights(orog_height, vert_grid):
    """
    Use the given orography and vertical discretization to compute the
    model heights.
    """
    nlon = len(orog_height.coord('longitude').points)
    nlat = len(orog_height.coord('latitude').points)

    # Skip the first level
    eta_theta = vert_grid._eta_theta[1:]
    first_const_rho_level = vert_grid._first_constant_rho

    nz = len(eta_theta)
    ztop = vert_grid._z_model_top
    
    height = numpy.zeros((nz, nlat, nlon))
    # Walk up to the first constant level
    for i in range(first_const_rho_level):
        height[i, :, :] = eta_theta[i] * ztop + orog_height.data * \
                numpy.square(
                        1 - eta_theta[i] / eta_theta[first_const_rho_level]
                        )

    # Remaining levels
    for i in range(first_const_rho_level, nz):
        height[i, :, :] = eta_theta[i] * ztop

    return height


def match_coord_system(source, target):
    """
    Set the coordinate systems on target to be the same as source, if not yet
    set.
    """
    target_dims = [coord.name() for coord in target.dim_coords]
    for coord in source.dim_coords:
        if (co_name := coord.name()) in target_dims:
            if target.coord(co_name).coord_system is not None:
                assert target.coord(co_name).coord_system == \
                        coord.coord_system, \
                        """Target coord already has a coordinate system that 
                        doesn't match the source coordinate system."""
            else:
                target.coord(co_name).coord_system = coord.coord_system

            
def vertical_interpolate(
        source,
        pressure_height,
        model_height,
        ):
    """
    Perform vertical interpolation of the source data to the target model
    heights, using the provided pressure heights.
    """

    nt_source, _, nlat, nlon = source.data.shape
    nt_press = pressure_height.shape[0]
    nz = model_height.shape[0]

    new_cube = iris.cube.Cube(
            numpy.zeros((nt_source, nz, nlat, nlon)),
            var_name=source.metadata.long_name,
            units=source.metadata.units,
            cell_methods=source.metadata.cell_methods,
            dim_coords_and_dims=[
                (source.coord('time'), 0),
                (iris.coords.DimCoord(
                    numpy.arange(1, nz + 1),
                    standard_name="model_level_number"
                    ), 1),
                (source.coord('latitude'), 2),
                (source.coord('longitude'), 3)
                ]
            )
    
    for (t, lat, lon) in numpy.ndindex(nt_source, nlat, nlon):
        local_pres_height = pressure_height.data[numpy.mod(t, nt_press), :, lat]
        local_source = source.data[t, :, lat, lon]
        local_model_height = model_height[:, lat, lon]

        interpolator = scipy.interpolate.interp1d(
                local_pres_height,
                local_source,
                kind='linear',
                bounds_error=False,
                fill_value=local_source[0]
                )

        new_cube.data[t, :, lat, lon] = interpolator(local_model_height)

    return new_cube


if __name__ == '__main__':
    args = parse_args()

    main(
            args.ozone_path,
            args.temperature_path,
            args.orography_file,
            args.vertical_disc,
            args.output,
            args.begin,
            args.end,
            args.zonal,
            args.climatological,
            args.extend,
            args.netcdf_only
            )
