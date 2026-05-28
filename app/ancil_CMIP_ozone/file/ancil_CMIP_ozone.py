#!/usr/bin/env python

#example usage...  python2.7 ozone_cmip6_ancillary_checked.py  -b 2050 -e 2099 -o /nesi/project/niwa00013//williamsjh/CMIP6/ozone-forcing/2019/trunk//OzoneConc/output/SSP5-8.5/ --levelfile /opt/niwa/um_sys/um/vn10.3/ctldata/vert/vertlevs_L85* --orogfile $UMDIR/ancil/atmos/n96e/orca025/orography/globe30/v7/qrparm.orog  --tempfile /nesi/project/niwa00013/williamsjh/CMIP6/ozone-forcing/2019/TASK-FORCE/zmta_input4MIPs_Temperature_CMIP_UReading-CCMI-SSP585-1-0_gn_201501-210012.nc --ozonefile /nesi/project/niwa00013/williamsjh/CMIP6/ozone-forcing/2019/jasmin/gws/nopw/j04/cmip6_prep_vol1/cmip6_ancils/data/input4MIPs_2018-12-18/input4MIPs/CMIP6/ScenarioMIP/UReading/UReading-CCMI-ssp585-1-0/atmos/mon/vmro3/gn/v20181101/vmro3_input4MIPs_ozone_ScenarioMIP_UReading-CCMI-ssp585-1-0_gn_205001-209912.nc --extend n --preind False --o3inputyr0 2050 --tinputyr0 2015
#
# get or find the CMIP6 ozone data, temperature on pressure etc
# manipulate this data to produce what is needed as input - subset of years, mean over certain years etc
# use the code to regrid this data to the required UM grid - horizontal and vertical
# write out the ancil file (as netcdf or ancil format?)
# need to know: years, timeseries, timeslice (or monthly climatology), extend time axis
#               UM model resolution (horiz + vert), orog file


'''
NAME:
    ozone_cmip6_ancillary_checked

DESCRIPTION:
    Generate the CMIP6 ozone ancillary for 1950-2014
    Based on idl code from Olaf Morgenstern
    Creates both 3D (tzyx) and 2D zonal mean (tzy1) ancillaries directly from netcdf to ancillary using ants

    Method:
        Set up 3D cube for new model grid based on resolution required
            use orography ancillary as standard for this model resolution (read from supercomputer)
        Read in zonal mean temperature-pressure file supplied by Michaela Hegglin
        Calculate pressure level boundaries
        Calculate heights of pressure levels - mid-point and edges

        Read in the CMIP6 ozone
        Rescale to MMR
        Read namelist file of eta_theta levels (scaled model levels)
        Calculate height of model levels scaled by orography

        Linearly interpolate ozone on pressure levels to model latitude/longitude
        Linearly interpolate height of pressure level edges to model latitudes
        Linearly interpolate each column of ozone using the pressure-height coordinates
            if the model level height goes below the lowest pressure level height, simply set the lowest
            model level ozone equal to the lowest pressure level ozone

    Input ozone climatology (pre-industrial 1850, historic 1850-2014)
        66 pressure levels, 1000 to 0.0001
        144 longitudes
        96 latitudes

    Input zonal mean temperature monthly climatology:
        CMIP6 temperature-pressure (monthly mean), on same grid as ozone
        (zonal mean supplied currently)

    Issues:
        Want to capture metadata:
            The input orography file used for the interpolation
        Converting from ppmv (SPARC) to mmr (model) - see http://www-nwp.metoffice.com/~magmp/MA_config/sparcoz_interpol.pro
        mass ozone Mozo = 48.0 (kg kMol-1), mass air Mair = 28.966 (kg kMol-1)
        ppmv2mmr = (1.0e-6 * Mozo / Mair)


AUTHOR:
    Based on idl code interpolate_height.pro by Olaf Morgenstern

    CMIP5 versions of the ozone code were:
    Based on code http://www-nwp.metoffice.com/~magmp/MA_config/sparcoz_interpol.pro to convert from
    CMIP5-based pressure level ozone (in ppmv)
    Based on code http://www-hc/~hadmn/VersionControl/utilities_088/v01/create_L85_ozone_ppfields.pro
    Malcolm Roberts (hadom)

    Expanded and generalised by Jonny Williams, NIWA, New Zealand

LAST MODIFIED:
    2021-08-27

'''

import sys
import matplotlib
import ants
import iris.coord_categorisation as icc
import os, glob
import numpy as np
import iris, cf_units, iris.analysis, iris.fileformats
import datetime, time
from scipy import interpolate
import matplotlib.pyplot as plt
import iris.plot as iplt

import argparse

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('filepaths', type=str, nargs='+',
                        help='Source data filepaths, relative to source-root.')

    _end_help = 'Compulsory end year for the processing.'
    parser.add_argument('-e', '--end', required=True, help=_end_help, type = int)

    _begin_help = 'Compulsory begin year for the processing.'
    parser.add_argument('-b', '--begin', required=True, help=_begin_help, type = int)

    _sources_help = 'sources'
    parser.add_argument('-s', '--sources', help=_sources_help, type = str)

    _output_help = 'Filename to write the result to or file path to write out multiple paths.'
    parser.add_argument('-o', '--output_dir', required=True,  help=_output_help, type = str)

    _levelfile_help = 'File describing the vertical levels in the atmosphere.'
    parser.add_argument('--levelfile', required=True,  help=_levelfile_help, type = str)

#    _ozonefile_help = 'File containing the input ozone data.'
#    parser.add_argument('--ozonefile', required=True,  help=_ozonefile_help, type = str, nargs='+')

    _target_grid_orog_help = 'File containing the orography data at the relevant resolution.'
    parser.add_argument('--target-grid-orog', required=True,  help=_target_grid_orog_help, type = str)

    _tempfile_help = 'File containing the input temperature data.'
    parser.add_argument('--tempfile', required=True,  help=_tempfile_help, type = str)
    
    #_preind_help = 'Flag determining whether this is a preindustrial (always year 0) run.'
    #parser.add_argument('--preind', required=True,  help=_preind_help, type = str)

    _timeaxis_help = 'Time axis of ancil - Timeseries, Timemean, Timeslice'
    parser.add_argument('--timeaxis', required=True,  help=_timeaxis_help, type = str)

    _extend_data_help = 'Flag determining whether to add data at start and end of timeseries'
    parser.add_argument('--extend_data', required=True,  help=_extend_data_help, type = str)

    #_atmos_resol_help = 'Atmosphere resolution.'
    #parser.add_argument('--atmos_resol', required=True,  help=_atmos_resol_help, type = str)

    #_ocean_resol_help = 'Ocean resolution.'
    #parser.add_argument('--ocean_resol', required=True,  help=_ocean_resol_help, type = str)

    _calendar_help = 'Model calendar.'
    parser.add_argument('--calendar', required=True,  help=_calendar_help, type = str)

    #_o3inputyr0_help = 'First year of data in the input ozone data'
    #parser.add_argument('--o3inputyr0', required=True,  help=_o3inputyr0_help, type = int)

    #_tinputyr0_help = 'First year of data in the input temperature data'
    #parser.add_argument('--tinputyr0', required=True,  help=_tinputyr0_help, type = int)

    #Help text partly copied from https://code.metoffice.gov.uk/trac/ancil/browser/contrib/trunk/AerosolChemistryEmissions/CMIP6/bin/preproc_emiss_cmip6.py

    #_extend_help = 'Extend the time series [y = yes, anything else for false]? This is required, for example, if you are running from 18500101 and your input data begins at the same time. The UM needs data before and after the times of interest in order to interpolate the data over time.'
    #parser.add_argument('--extend', required=True,  help=_extend_help, type = str)

    args = parser.parse_args()

    return args

def zonal_callback(cube, field, filename):
    '''
    Currently the way to include information needed to make (zonal mean) ancillary
    '''
    for nc, coord in enumerate(cube.dim_coords):
        if 'longitude' in coord.name():
            coord_longitude = nc

    if cube.shape[-1] == 1:
        print('zonal mean input file ',cube)

        xc = iris.coords.DimCoord([180, ],
                                  bounds=[[0, 360], ],
                                  standard_name='longitude',
                                  units='degrees')
        cube.remove_coord('longitude')
        cube.add_dim_coord(xc, coord_longitude)
        cm = iris.coords.CellMethod(method='mean', coords='longitude')
        cube.add_cell_method(cm)

    try:
        cube.remove_coord('model_level_number')
    except:
        cube.remove_coord('level_height')

    #cube.remove_coord('hybrid_ht')
    # 3 hardcoded here - be sure it's the right coordinate

    # Cell methods required for populating LBPROC for each field
    cm = iris.coords.CellMethod(method='mean', coords='time')
    cube.add_cell_method(cm)
    try:
        ants.utils.guess_bounds(cube.coord('time'))
    except:
        pass

def setup_model_grid(level_namelist_file, orog_height, year, calendar):
    '''
    return a iris cube with:
        the meridional number of grid points based on the half-grid of the zonal component
        (N in UM speak)
    '''
    latitudes = orog_height.coord('latitude').points
    longitudes = orog_height.coord('longitude').points

    vertical_levels = ants.load_grid(level_namelist_file)

    level_height = vertical_levels.coord('level_height').points
    #print('level_height ',level_height, len(level_height))
    model_levs = np.arange(len(level_height))+1

    xcoord = orog_height.coord('longitude')
    ycoord = orog_height.coord('latitude')
    zcoord = iris.coords.DimCoord(level_height, long_name = zcoord_name, units='m')

    # Create a new time coordinate
    time_pts = []
    if calendar == '360_day':
        calendar_unit = cf_units.CALENDAR_360_DAY
    else:
        raise Exception('Calendar not 360_day'+calendar)

    n_months = 12
    data = np.zeros((n_months, len(model_levs), len(latitudes), len(longitudes)))
    for month in range(1,13):
        i = datetime.datetime(year, month, 15,0,0,0)
        time_points = cf_units.date2num(i, 'hours since 1850-01-01 00:00:00',calendar=calendar)
        time_pts.append(time_points)
    time_coord = iris.coords.DimCoord(time_pts, standard_name='time',
                     units=cf_units.Unit('hours since 1850-01-01 00:00:00',
                                 calendar=calendar_unit))
    cube=iris.cube.Cube(data, var_name = 'Ozone', long_name='SPARC OZONE', units = '')
    cube.add_dim_coord(zcoord,1)
    cube.add_dim_coord(ycoord,2)
    cube.add_dim_coord(xcoord,3)
    cube.add_dim_coord(time_coord,0)

    return cube

def vertical_resolution_check(zlevels):
    if zlevels != 38 or zlevels != 85:
        sys.exit('unrecognised vertical levels '+zlevels)
    else:
        pass

def read_cmip6_ozone(filename, o3input_startyear):
    '''
    Read ozone data into a cube
    For the initial CMIP6 pre-industrial file, need to set up a proper calendar time coordinate
    (assume 1850 here), since iris does not like month units
    Read the ozone data in the directory+filename path
    Use year to extract the appropriate year if needed

    Arguments:
        directory + filename - path to ozone dataset
        data_year - integer - the initial year in the dataset (to reset the calendar) - default 1850
        extract_years - list of 2 integers - years to extract from dataset (inclusive)

    '''
    data_year = o3input_startyear
    cube_ozone_list = iris.load(filename)
    print(cube_ozone_list)
    from iris.experimental.equalise_cubes import equalise_attributes
    equalise_attributes(cube_ozone_list)
    cube_ozone = cube_ozone_list.concatenate_cube()
    print('cube_ozone ',cube_ozone.shape)
    print('cube_ozone time ',cube_ozone.coord('time'))

    nyears = int(cube_ozone.shape[0] / 12)
    years = list(range(data_year, data_year+nyears,1))
    # fix the time coordinate in the ozone file
    time_pts = []
    for year in years:
        for month in range(1,13):
            i = datetime.datetime(year, month, 15,0,0,0)
            tunit = 'hours since '+str(data_year)+'-01-01 00:00:00'
            time_points = cf_units.date2num(i, tunit, calendar='gregorian')
            time_pts.append(time_points)
    time_coord = iris.coords.DimCoord(time_pts, standard_name='time',
                     units=cf_units.Unit(tunit,
                                 calendar=cf_units.CALENDAR_GREGORIAN))
    cube_ozone.remove_coord('time')
    cube_ozone.add_dim_coord(time_coord, 0)

    print(cube_ozone.coord('time'))

    icc.add_year(cube_ozone, 'time', 'year')
    icc.add_month_number(cube_ozone, 'time', 'month')
    o3input_endyear = cube_ozone.coord('year').points[-1]

    return cube_ozone, o3input_endyear

def read_cmip6_temperature(filename, tinputyr0):
    '''
    read in the cmip6 supplied pre-industrial, zonal mean temperature file in netcdf
    For the initial CMIP6 pre-industrial file, need to set up a proper calendar time coordinate
    (assume 1850 here) - iris does not like months
    Read the ozone data in the directory+filename path
    Use year to extract the appropriate year if needed
    Arguments:
        directory + filename - path to ozone dataset
        data_year - integer - the initial year in the dataset (to reset the calendar) - default 1850
        extract_years - list of 2 integers - years to extract from dataset (inclusive)
    '''
    data_year = tinputyr0
    cube_tmzm = iris.load_cube(filename)
    nyears = int(cube_tmzm.shape[0] / 12)
    print('nyears ',nyears, cube_tmzm.shape[0])
    years = list(range(data_year, data_year+nyears,1))

    # fix the time coordinate in the temperature file
    time_pts = []
    for year in years:
        for month in range(1,13):
            i = datetime.datetime(year, month, 15,0,0,0)
            tunit = 'hours since '+str(data_year)+'-01-01 00:00:00'
            time_points = cf_units.date2num(i, tunit,calendar='gregorian')
            time_pts.append(time_points)
    time_coord = iris.coords.DimCoord(time_pts, standard_name='time',
                     units=cf_units.Unit(tunit,
                                 calendar=cf_units.CALENDAR_GREGORIAN))
    cube_tmzm.remove_coord('time')
    cube_tmzm.add_dim_coord(time_coord, 0)

    icc.add_year(cube_tmzm, 'time', 'year')
    icc.add_month_number(cube_tmzm, 'time', 'month')
    return cube_tmzm

def zonal_mean_ozone(cube):
    cube_zonal = cube.collapsed('longitude', iris.analysis.MEAN)
    cube_zonal_lon = iris.util.new_axis(
            cube_zonal,
            scalar_coord=iris.coords.DimCoord(
                [360.0],
                standard_name='longitude',
                circular=True,
                units='degrees',
                coord_system=iris.coord_system.GeogCS(
                    iris.fileformats.pp.EARTH_RADIUS
                    )
                )
            )
    )

    cube_zonal_lon.transpose(cube_zonal.shape, (1,))

    return cube_zonal_lon

def interpolate_to_new_latitudes(from_cube, to_cube):
    interpolator = iris.analysis.Linear(extrapolation_mode='linear').interpolator(from_cube, ['latitude'])
    # these are the new latitudes to interpolate to
    new_latitude = to_cube.coord('latitude').points
    # now interpolate to create new cube
    result = interpolator([new_latitude])
    return result

def read_in_model_levels(level_namelist_file):
    '''
    Read in the model levels from a dataset name containing the inlevs string
    '''
    vertical_levels = ants.load_grid(level_namelist_file)

    return vertical_levels.coord('level_height').points

def interpolate_column(from_data_cube, z2_pressure_height, ht_model_height_coord, to_cube):
    '''
    interpolate data values held in 'from_data_cube' from its pressure levels to model levels

    Inputs:
        from_data_cube: pressure level ozone data, nt, npressure, latitude, longitude
        z2_pressure_height: height of pressure levels
        ht_model_height_coord: height of model hybrid coordinate
        to_cube: reference cube on output grid
    '''
    print('interp ',to_cube)
    #new_cube = to_cube.copy()
    new_cube = iris.cube.CubeList()

    for imon in range(1,13):
        print('month ',imon)
        #month_constraint = iris.Constraint(coord_values = {'month_number' : lambda l : l == imon})
        #sub_to_cube = to_cube.extract(month_constraint)
        #sub_from_cube = from_data_cube.extract(month_constraint)
        sub_to_cube = to_cube[imon-1]
        sub_from_cube = from_data_cube[imon-1]
        cube_interp = sub_to_cube.copy()
        for j, lat in enumerate(sub_to_cube.coord('latitude').points):
            press_height = z2_pressure_height.data[imon-1, :, j]
            for i, longit in enumerate(sub_to_cube.coord('longitude').points):
                model_height = ht_model_height_coord[:, j, i]
                model_data = sub_from_cube.data[:, j, i]

                # use bottom value for any points where new levels lie outside the range of the old levels
                # (i.e. would need extrapolation) - without this part, this version of interp1d simply
                # fails (more recent versions have options to better cope)
                fill_value = sub_from_cube.data[0, j, i]
                interpfunc = interpolate.interp1d(press_height, model_data, kind='linear', \
                                                  bounds_error=False, fill_value=fill_value)
                result=interpfunc(model_height)

                #new_cube.data[imon, :, j, i] = result
                cube_interp.data[:, j, i] = result
        cube_interp = iris.util.new_axis(cube_interp, 'time')
        new_cube.append(cube_interp)
    new_cube = new_cube.concatenate_cube()
    new_cube.var_name = 'Ozone'
    new_cube.long_name = 'Ozone MMR'
    return new_cube

def make_column_cube(data, coord, coord_name):
    '''
    Construct a cube from the data array
    '''
    zcoord=iris.coords.DimCoord(coord, long_name=coord_name, units='')
    cube=iris.cube.Cube(data, var_name = 'Ozone', long_name='SPARC OZONE', units = '')
    cube.add_dim_coord(zcoord,0)
    return cube

def add_pp_information(cube):
    '''
    Add pp specific header information so that ancillary creation can work
    cube - input cube of data
    variable - variable name (from which Met Office PP information will be derived)
    Returns an iterable
    '''
    cube.attributes['STASH'] = iris.fileformats.pp.STASH.from_msi('m01s00i060')
    for cube, field in iris.fileformats.pp.as_pairs(cube):
        # post process the PP field, prior to saving
        field.lbfc = 453
        # This overrides the existing setting - is the new value different though?
        field.lblev = function_height_level(cube.coord(zcoord_name).points)
#        cube.standard_name = 'surface_temperature'
        yield field

def function_height_level(height):
    level_height = read_in_model_levels(INLEV)
    miss = np.where(height == level_height)
    print(miss[0]+1)
    return miss[0]+1

def derive_pressure_heights(cube, plevs):
    '''
    Derive the height on pressure level middle and edges
    cube - iris cube of (zonal mean) temperature on pressure levels
    plevs - pressure levels 

    Output:
        cube of pressure heights at edges
    '''
    dims = cube.shape
    ntime = dims[0]; npres = dims[1]; nlat = dims[2]

    integrant = R_gas * cube / Grav

    # calculate pressure boundary values - assume bottom is sea level at 1013 and top is 1e-5
    pboundary = np.zeros(npres+1)
    pboundary[0] = Sealevpr
    pboundary[-1] = 1.0e-5
    pboundary[1:-1] = np.sqrt(plevs[0:-1]*plevs[1:])

    log_pboundary = np.log(pboundary)
    dlogp = np.log(pboundary[0:-1] / pboundary[1:])
    frac = (log_pboundary[0:-1] - np.log(plevs)) / (log_pboundary[0:-1] - log_pboundary[1:])

    # z = height at edge of box
    # z1 = height at centre of box
    z = np.zeros((ntime, npres+1, nlat))
    z1 = cube.copy()

    # integrate upwards
    for i in range(npres):
        z[:, i+1, : ] = z[:, i, :] + dlogp[i] * integrant.data[:, i, :]
        z1.data[:, i, : ] = z[:, i, :] + frac[i] * dlogp[i] * integrant.data[:, i, :]

    return z1

def compare_temp_ozone_years(temp_on_press, ozone_on_press, override = False):
    '''
    Check that years in ozone and temperature files are the same
    '''
    years_ozone = ozone_on_press.coord('year').points
    years_temp = temp_on_press.coord('year').points

    if len(years_ozone) != len(years_temp):
        raise Exception('Looks like years in ozone and temperature not the same')

    for yo, yt in zip(years_ozone, years_temp):
        if yo != yt:
            if override:
                print('WARNING, years in temp and ozone do not agree ')
            else:
                raise Exception('Looks like years in ozone and temperature not the same')

def calculate_model_heights(orog_height, level_namelist_file):
    '''
    Calculate model hybrid heights

    Output:
        Numpy array of (model level, lat, long) of heights scaled by orography
    '''
    dims_orog = orog_height.shape
    nlonorog = dims_orog[1]; nlatorog = dims_orog[0]

    # file containing model level information (eta-theta and eta_rho, model top etc)
    # use some ants code to get the information needed
    #namelist = ants.fileformats.namelist._read_namelist(level_namelist_file)['vertlevs']
    namelist = ants.fileformats.namelist._read_namelist(level_namelist_file)
    #print('namelist ',namelist)
    vertical_levels = ants.fileformats.namelist.VerticalLevels(namelist)
    # skip surface level
    eta_theta = vertical_levels._eta_theta[1:]
    first_constant_r_rho_level = vertical_levels._first_constant_rho

    nz = len(eta_theta)
    ztop = vertical_levels._z_model_top
    print ('eta_levels ',eta_theta, ztop)

    # height of model levels scaled by orography
    ht = np.zeros((nz, nlatorog, nlonorog))
    for i in range(0, first_constant_r_rho_level):
        ht[i, :, :] = eta_theta[i]*ztop + orog_height.data[:, :] * np.square(1.0 - eta_theta[i] / eta_theta[first_constant_r_rho_level])
    for i in range(first_constant_r_rho_level, nz):
        ht[i, :, :] = eta_theta[i]*ztop

    return ht

def read_cmip6_data(temp_file, ozone_file, massmix, tinputyr0, o3inputyr0, override=False):
    '''
    Read CMIP6 input data, ozone on pressure levels and temperature on pressure
    Convert ozone to mass mixing ratio if required
    Check that dates match between ozone and temperature
    '''
    temp_on_press_cmip6 = read_cmip6_temperature(temp_file, tinputyr0)
    ozone_on_press_cmip6, o3input_endyear = read_cmip6_ozone(ozone_file, o3inputyr0)
    print('ppmv2mmr ',ppmv2mmr)
    print('temp_on_press_cmip6 ',temp_on_press_cmip6)
    print('ozone_on_press_cmip6 ',ozone_on_press_cmip6)
    if massmix: ozone_on_press_cmip6 *= ppmv2mmr

    # check years are the same at this point
    compare_temp_ozone_years(temp_on_press_cmip6, ozone_on_press_cmip6, override=True)

    return temp_on_press_cmip6, ozone_on_press_cmip6, o3input_endyear

def work(temp_on_press_cmip6, ozone_on_press_cmip6, calendar, level_namelist_file, final_3dfile, final_2dfile):
    '''
    Process ozone data on pressure levels and produce it on model levels
    if we produce a timeseries, then loop through years, extract each year of data and interpolate
    if we produce a (time-mean) climatology, then just need the interpolation for one year
    '''
    # Read orography
    orog_height = iris.load_cube(args.target_grid_orog,'surface_altitude')

    # calculate pressure height cube
    plevs = temp_on_press_cmip6.coord('air_pressure').points
    z1_pressure_height = derive_pressure_heights(temp_on_press_cmip6, plevs)

    # interpolate z1_pressure_height to model latitudes
    z2_pressure_height = interpolate_to_new_latitudes(z1_pressure_height, orog_height)

    # calculate model heights
    ht_model_levels = calculate_model_heights(orog_height, level_namelist_file)

    # make sure coord_system is same else this does not work
    for coord in ['latitude','longitude']:
        ozone_on_press_cmip6.coord(coord).coord_system = orog_height.coord(coord).coord_system

    # set up cube lists to hold the year by year data
    full_ozone_cube_list = iris.cube.CubeList()
    full_ozone_zonal_cube_list = iris.cube.CubeList()

    year_cmip_begin = int(temp_on_press_cmip6.coord('year').points[0])
    year_cmip_end = int(temp_on_press_cmip6.coord('year').points[-1]+1)
    print('year_cmip_begin , year_cmip_end ',year_cmip_begin , year_cmip_end)

    # here the year range needs to be based on how many years (i.e. timeseries, climatology etc)
    for year in range(year_cmip_begin, year_cmip_end):
        file_ozone_3D_year = final_3dfile+'_yr'+str(year)+'.nc'
        file_ozone_zonal_year = final_2dfile+'_yr'+str(year)+'.nc'
        if not os.path.exists(file_ozone_zonal_year):
            print('year ',year)
            model_cube = setup_model_grid(level_namelist_file, orog_height, year, calendar)
            # interpolate ozone onto model grid
            year_constraint = iris.Constraint(coord_values = {'year' : lambda l : l == year})
            temp_year = temp_on_press_cmip6.extract(year_constraint)
            ozone_year = ozone_on_press_cmip6.extract(year_constraint)

            ozone_lat_long = ozone_year.regrid(orog_height, iris.analysis.Linear())
            z2_pressure_height_year = z2_pressure_height.extract(year_constraint)
            # interpolate ozone to model levels
            # 66 pressure levels
            # use ht_model_levels model heights
            ozone_model_height = interpolate_column(ozone_lat_long, z2_pressure_height_year, ht_model_levels, model_cube)

            ozone_model_height.standard_name = 'mass_fraction_of_ozone_in_air'
            # add a stashcode for conversion to ancillary
            ozone_model_height.attributes['STASH'] = iris.fileformats.pp.STASH.from_msi('m01s00i060')

            add_pp_information(ozone_model_height)
            ozone_zonal = zonal_mean_ozone(ozone_model_height)

            iris.save(ozone_model_height, file_ozone_3D_year, unlimited_dimensions = ['time'], netcdf_format='NETCDF4', zlib=True, complevel=2)
            iris.save(ozone_zonal, file_ozone_zonal_year, unlimited_dimensions = ['time'], netcdf_format='NETCDF4', zlib=True, complevel=2)
        
        ozone_model_height = iris.load_cube(file_ozone_3D_year)
        ozone_zonal = iris.load_cube(file_ozone_zonal_year)

        full_ozone_cube_list.append(ozone_model_height)
        full_ozone_zonal_cube_list.append(ozone_zonal)

    print('full_ozone_zonal_cube_list ',full_ozone_zonal_cube_list)
    full_ozone_cube = full_ozone_cube_list.concatenate_cube()

    full_ozone_zonal_cube = full_ozone_zonal_cube_list.concatenate_cube()
    print('full_ozone_zonal_cube ',full_ozone_zonal_cube)

    # save both files
    iris.save(full_ozone_cube, final_3dfile, unlimited_dimensions = ['time'], netcdf_format='NETCDF4', zlib=True, complevel=2)
    iris.save(full_ozone_zonal_cube, final_2dfile, unlimited_dimensions = ['time'], netcdf_format='NETCDF4', zlib=True, complevel=2)

def produce_ancillary_with_ants(final_3d_source, final_2d_source, start_year, end_year, timeaxis, extend_data, o3input_startyear, o3input_endyear, level_namelist_file):

    for fname in [final_2d_source, final_3d_source]:
        cube = ants.load_cube(fname, callback = zonal_callback)
        # Vertical coords from namelist
        vertical_levels = ants.load_grid(level_namelist_file)
        cube.add_dim_coord(vertical_levels.coord('model_level_number'), 1)
        cube.add_aux_coord(vertical_levels.coord('sigma'), 1)
        cube.add_aux_coord(vertical_levels.coord('level_height'), 1)

        try:
            cube.coord('time').guess_bounds()
        except:
            pass

        try:
            cubes.remove_coord('year')
            cubes.remove_coord('month')
        except:
            pass

        #Following `if statement` copied from partly copied from https://code.metoffice.gov.uk/trac/ancil/browser/contrib/trunk/AerosolChemistryEmissions/CMIP6/bin/preproc_emiss_cmip6.py
        print('extend_data, start year ',extend_data, start_year, o3input_startyear)
        if extend_data:
            tc = cube.coord(axis='t')
            if end_year == o3input_endyear:
                # Duplicate the final year data if end is the year after the raw
                # source data runs out.  UM needs data before and after the times of
                # interest in order to interpolate the data over time.
            
                length_one_year = tc.points[-1] - tc.points[-13]
                extra_year = cube[-12:].copy()
                extra_year_tc = extra_year.coord(axis='t')
                extra_year_tc.points = extra_year_tc.points + length_one_year
                extra_year_tc.bounds = extra_year_tc.bounds + length_one_year
                cubelist = iris.cube.CubeList((cube, extra_year))
                cube = cubelist.concatenate_cube()

            if start_year == o3input_startyear:
                # And the same for the start year
                length_one_year = tc.points[12] - tc.points[0]
                extra_year = cube[:12].copy()
                extra_year_tc = extra_year.coord(axis='t')
                extra_year_tc.points = extra_year_tc.points - length_one_year
                extra_year_tc.bounds = extra_year_tc.bounds - length_one_year
                cubelist = iris.cube.CubeList((extra_year, cube))
                cube = cubelist.concatenate_cube()

        elif not extend_data: 
            pass

        if timeaxis == 'Timemean':
            if cube.coord('time').has_bounds():
                cube.coord('time').bounds = None

            cube.coord('time').guess_bounds()

            tc = cube.coord('time')
            tc.points = np.around(tc.points)
            tc.bounds = np.around(tc.bounds)

            # do we need to set the year to 0000 (I guess we can't do that with datetime)?
            #tunit = cube.coord('time').units
            #new_date = []
            #for imon in range(1,13):
            #    new_date.append(tunit.date2num(datetime.datetime(1950,imon,16,0,0,0)))
            #    cube.coord('time').points = new_date


        outfile_anc = fname[:-3]+'_ants.anc'
        ants.save(cube, outfile_anc, saver='ancil')

def fix_calendar_1950(cube):
    # For the first and last time points, the calendar needs adjusting to be the previous month
    # need to set the calendar to yyyymmm_out
    cube.coord('time').bounds = None
    tunit = cube.coord('time').units
    new_date = []
    for imon in range(1,13):
        new_date.append(tunit.date2num(datetime.datetime(1950,imon,16,0,0,0)))
    cube.coord('time').points = new_date

def produce_cyclic_ancillary_with_ants(source_file):
    cubes = ants.load_cube(source_file, callback = zonal_callback)  # Source filename
    try:
        cubes.remove_coord('year')
        cubes.remove_coord('month')
    except:
        pass
    if cubes.coord('time').has_bounds():
        cubes.coord('time').bounds = None

    cubes.coord('time').guess_bounds()

    tc = cubes.coord('time')
    tc.points = np.around(tc.points)
    tc.bounds = np.around(tc.bounds)
    ants.save(cubes, source_file[:-3]+'_ants.anc', saver='ancil')  # Output filename

def define_constants(match_olaf=True):
    '''
    Define constants used in code
    Parameters:
      match_olaf: logical: Match the values used in original CMIP6 code
    '''

    global MASS_OZONE, MASS_AIR, Navo, Ride, R_gas, Sealevpr, ppmv2mmr, Grav

    MASS_OZONE = 48.0 # (kg kMol^-1)
    MASS_AIR = 28.97 # (kg kMol^-1) as in the UM 6.1, see M_AIR
    Navo=6.022*1e4 # (10^22 molecules kMol^-1)
    Ride=8314.        # (J kg^-1 K^-1)
    R_gas = 287.058

    # use constants to match Olaf's code
    
    if match_olaf:
        Sealevpr=1013.  # (hPa) and is assumed constant at all latitudes
        ppmv2mmr = 1.65
        Grav=9.81      # (m s^-2) same value as for the UM 6.1
    else:
        Sealevpr=1013.25  # (hPa) and is assumed constant at all latitudes
        ppmv2mmr=(MASS_OZONE/MASS_AIR)
        Grav=9.80665      # (m s^-2) same value as for the UM 6.1

    MASSMIX = True  # convert to mass mixing ration (from ppmv)

def preprocess_input_data(cmip6_temp, cmip6_ozone, year_start, year_end, extend_data, timeaxis, massmix, tinput_startyear):
    '''
    Preprocess inputs, e.g. extract time period, take average over time period, etc
    We subset the time period here. If we need to extend the data, then we can add to the time period.
    However, if we are at the edges of the data, then do we need to extend in a different way?
    We also take the mean if climatology is true
    :param str indir
    :param str outdir
    '''
    # Read CMIP6 input data
    temp_on_press_cmip6, ozone_on_press_cmip6, o3input_endyear = read_cmip6_data(cmip6_temp, cmip6_ozone, massmix, tinput_startyear, o3input_startyear, override=True)

    if timeaxis == 'Timeseries':
        # need to extract the data for the times given, but may need to extend
        if extend_data:
            year_start_extract = year_start - 1
            year_end_extract = year_end + 1
        else:
            year_start_extract = year_start
            year_end_extract = year_end
    else:
        year_start_extract = year_start
        year_end_extract = year_end

    year_constraint = iris.Constraint(coord_values = {'year' : lambda l : year_start_extract <= l <= year_end_extract})
    temp = temp_on_press_cmip6.extract(year_constraint)
    temp_on_press_cmip6 = temp
    ozone = ozone_on_press_cmip6.extract(year_constraint)
    ozone_on_press_cmip6 = ozone

    if timeaxis == 'Timemean':
        # want monthly means over whole period
        temp = temp_on_press_cmip6.aggregated_by('month', iris.analysis.MEAN)
        temp_on_press_cmip6 = temp
        ozone = ozone_on_press_cmip6.aggregated_by('month', iris.analysis.MEAN)
        ozone_on_press_cmip6 = ozone

    return temp_on_press_cmip6, ozone_on_press_cmip6, o3input_endyear

def set_output_filenames(year_start, year_end, output_dir, timeaxis):
    if timeaxis == 'Timemean':
        OUTPUT_3D_FILE = 'mmro3_monthly_CMIP6v2_clim_3d_{}.nc'
        OUTPUT_2D_FILE = 'mmro3_monthly_CMIP6v2_clim_zonalmn_{}.nc'
        NYEARS = 1
    elif timeaxis == 'Timeseries' or timeaxis == 'Timeslice':
        OUTPUT_3D_FILE = 'mmro3_monthly_CMIP6v2_3d_{}.nc'
        OUTPUT_2D_FILE = 'mmro3_monthly_CMIP6v2_zonalmn_{}.nc'

    year_period = str(year_start)+'_'+str(year_end)

    # calculate the model ozone from the CMIP6 input on pressure levels and
    # temperature on pressure levels
    final_3dfile = os.path.join(output_dir, OUTPUT_3D_FILE.format(year_period))
    final_2dfile = os.path.join(output_dir, OUTPUT_2D_FILE.format(year_period))
    print(final_2dfile, final_3dfile)

    return final_3dfile, final_2dfile

def set_logical(value):
    if isinstance(value, str):
        if 'Y' in value or 'y' in value or 'True' in value or 'true' in value:
            return True
        else:
            return False
    else:
        return value

def check_timeaxis(value):
    if value == 'Timeseries' or value == 'Timemean' or value == 'Timeslice':
        return value
    else:
        raise Exception('Timeaxis value not one of Timeseries, Timemean, Timeslice')

if __name__ == '__main__':
    print (sys.argv[1:])
    args = parse_args()

    # horizontal resolution
    #atmos_resol = args.atmos_resol
    #ocean_resol = args.ocean_resol

    #cmip6_ozone_files = args.ozonefile
    cmip6_ozone_files = args.filepaths
    print(cmip6_ozone_files)
    cmip6_temp_file = args.tempfile

    year_start = args.begin
    year_end = args.end
    if year_end < year_start:
        raise Exception('year_end must be after year_start ')
    calendar = args.calendar
    
    # Need to consider timeseries, climatology, timeslice (same as timeseries?)
    # need a default for the timeseries e.g. 1850-2014, or has to be input
    # How to specify number of years to use for mean, or just assume n
    print('timeaxis ',args.timeaxis, args.extend_data)
    timeaxis = check_timeaxis(args.timeaxis)
    extend_data = set_logical(args.extend_data)
    print('timeseries ',timeaxis, type(timeaxis), extend_data)
    if extend_data and timeaxis != 'Timeseries':
        raise Exception('Cannot extend the data if not Timeseries')
    if timeaxis == 'Timeslice' and year_start != year_end:
        raise Exception('For Timeslice data, start and end year are the same')
                        
    level_file = args.levelfile
    print(level_file)

    o3input_startyear = 1850
    tinput_startyear = 1850
    massmix = True

    zcoord_name = 'model_level_number'

    output_dir = args.output_dir
    if output_dir is None:
        raise exception('Output_dir is not defined')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    finaldir = output_dir

    # define the constants to be used in the code
    define_constants()

    # set the (temporary) output filenames
    final_3dfile, final_2dfile = set_output_filenames(year_start, year_end, output_dir, timeaxis)

    # preprocess the input data
    temp_on_press_cmip6, ozone_on_press_cmip6, o3input_endyear = preprocess_input_data(cmip6_temp_file, cmip6_ozone_files, year_start, year_end, extend_data, timeaxis, massmix, tinput_startyear)
    print('temp_on_press ',temp_on_press_cmip6)
    print('ozone_on_press ',ozone_on_press_cmip6)

    # do the interpolation
    work(temp_on_press_cmip6, ozone_on_press_cmip6, calendar, level_file, final_3dfile, final_2dfile)

    # produce the final ancillary files
    produce_ancillary_with_ants(final_3dfile, final_2dfile, year_start, year_end, timeaxis, extend_data, o3input_startyear, o3input_endyear, level_file)
