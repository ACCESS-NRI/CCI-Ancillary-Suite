# Convert the specified ocean mesh into a land mask and land fractions
import re
import xarray
import esmpy
import numpy
import argparse
import ants

def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--ocean-mesh-file',
                        required=True,
                        help='ACCESS-OM3 mesh to process.'
                        )

    parser.add_argument('--output',
                        required=True,
                        help='Path to write the mask and land fractions to.'
                        )

    parser.add_argument('--atm-resolution',
                        required=True,
                        help='String describing the atmosphere resolution e.g. n96e, n512.'
                        )

    return parser.parse_args()

def parse_resolution_string(atm_res):
    """Parse the string given for the atmosphere resolution and return the
    number of longitude/latitude points."""

    if m := re.match('n(\d+)e$', atm_res):
        grid_type = 'EG'
    elif m := re.match('n(\d+)e$', atm_res):
        grid_type = 'ND'
    else:
        raise Exception(f'String {atm_res} is not a valid UM atmosphere resolution.')

    res = int(m.group(1))
    assert res % 2 == 0, 'Resolution must be a multiple of 2.'
    nlon = 2 * res
    nlat = 3 * res // 2

    return nlon, nlat

def load_ocean_data(ocn_mesh_fp):
    mesh_ds = xarray.load_dataset(ocn_mesh_fp, engine='netcdf4')

    ocn_mask = mesh_ds.elementMask.data

    ocn_mesh = esmpy.api.mesh.Mesh(
    filename=ocn_mesh_fp,
    filetype=esmpy.api.constants.FileFormat.ESMFMESH
    )

    return ocn_mesh, ocn_mask

def create_um_mesh(nlat, nlon):
    # make UM mesh file from the required latitudes and longitudes

    num_elements = nlon * nlat

    element_ids = numpy.arange(1, num_elements + 1)
    element_types = numpy.array([
      esmpy.MeshElemType.TRI] * nlon + [esmpy.MeshElemType.QUAD] * ((nlat - 2) * nlon) + [esmpy.MeshElemType.TRI] * nlon
                          )
    element_lon = numpy.zeros(num_elements)
    element_lon[:] = ((element_ids - 1) % nlon + 0.5) * (360 / nlon)

    element_lat = numpy.zeros(num_elements)
    element_lat[:] = ((element_ids - 1) // nlon + 0.5) * (180 / nlat)  - 90

    element_coords = numpy.zeros((num_elements, 2))
    element_coords[:, 0] = element_lon
    element_coords[:, 1] = element_lat

    dx = 360 / nlon
    dy = 180 / nlat
    pi_over_180 = numpy.pi / 180
    element_areas = dx * pi_over_180 * (
      numpy.sin((element_lat + 0.5 * dy) * pi_over_180) - numpy.sin((element_lat - 0.5 * dy) * pi_over_180)
    )

    num_nodes = nlon * (nlat - 1) + 2
    node_ids = numpy.arange(1, num_nodes + 1)

    node_lon = numpy.zeros(num_nodes)
    node_lon[0] = 0.0
    node_lon[-1] = 0.0
    node_lon[1:-1] = element_lon[:-nlon] - 0.5 * (360 / nlon)

    node_lat = numpy.zeros(num_nodes)
    node_lat[0] = -90.0
    node_lat[-1] = 90.0
    node_lat[1:-1] = element_lat[:-nlon] + 0.5 * (180 / nlat)

    node_coords = numpy.zeros((num_nodes, 2))
    node_coords[:, 0] = node_lon
    node_coords[:, 1] = node_lat

    element_conn = []

    iy = 0
    for ix in range(nlon):
      south_pole = 1
      north_west = iy * nlon + 2 + ix
      north_east = iy * nlon + 2 + (ix + 1) % nlon
      conn = [north_west, south_pole, north_east]

      element_conn.extend(conn)
      # print(conn)


    for iy in range(1, nlat - 1):
      for ix in range(nlon):
          north_west = iy * nlon + 2 + ix
          north_east = iy * nlon + 2 + (ix + 1) % nlon

          south_west = (iy - 1) * nlon + 2 + ix
          south_east = (iy - 1) * nlon + 2 + (ix + 1) % nlon

          conn = [north_west, south_west, south_east, north_east]
          element_conn.extend(conn)

    iy = nlat - 1
    for ix in range(nlon):
      north_pole = num_nodes
      south_west = (iy - 1) * nlon + 2 + ix
      south_east = (iy - 1) * nlon + 2 + (ix + 1) % nlon
      conn = [north_pole, south_west, south_east]

      element_conn.extend(conn)

    element_conn = numpy.array(element_conn)

    nodeCoords = xarray.DataArray(
      dims=("nodeCount", "coordDim"),
      data=node_coords,
      attrs=dict(units="degrees")
    )

    elementConn = xarray.DataArray(
      dims=("connectionCount"),
      data=element_conn.astype(numpy.int32),
      attrs=dict(long_name="Node indices that define the element connectivity")
    )

    numElementConn = xarray.DataArray(
      dims=("elementCount"),
      data=element_types.astype(numpy.int32),
      attrs=dict(long_name="Number of nodes per element")
    )

    centerCoords = xarray.DataArray(
      dims=("elementCount", "coordDim"),
      data=element_coords,
      attrs=dict(units="degrees")
    )

    elementArea = xarray.DataArray(
      dims=("elementCount"),
      data=element_areas,
    )

    um_mesh_ds = xarray.Dataset(
      dict(nodeCoords=nodeCoords, elementConn=elementConn, numElementConn=numElementConn, centerCoords=centerCoords, elementArea=elementArea),
      attrs=dict(gridType="unstructured mesh")
    )

    um_mesh = esmpy.api.mesh.Mesh(
    filename='um_mesh.nc',
    filetype=esmpy.api.constants.FileFormat.ESMFMESH
    )

    return um_mesh

def create_landfracs(um_mesh, ocn_mesh, ocn_mask):
    # Create the land fracton array from the ocean mesh

    atm_field = esmpy.Field(um_mesh, meshloc=esmpy.api.constants.MeshLoc.ELEMENT)
    atm_field.data[:] = 0.0


    ocn_field = esmpy.Field(ocn_mesh, meshloc=esmpy.api.constants.MeshLoc.ELEMENT)
    ocn_field.data[:] = ocn_mask

    ocn_to_atm_cons = esmpy.api.regrid.Regrid(
      ocn_field, atm_field,
      unmapped_action=esmpy.api.constants.UnmappedAction.IGNORE,
      regrid_method=esmpy.api.constants.RegridMethod.CONSERVE,
      norm_type=esmpy.api.constants.NormType.DSTAREA, factors=True
    )

    new_ocn_frac = atm_field.data.reshape((nlat, nlon))

    new_land_frac = 1.0 - new_ocn_frac

    new_land_frac[new_land_frac < 0.01] = 0.0
    new_land_frac[new_land_frac > 1.0] = 1.0

    return new_land_frac


def save_landfracs(nlat, nlon, land_frac, out_fp):
    # Write the land fractions to NetCDF and as UM ancillary file
    dx = 360 / nlon
    dy = 180 / nlat

    lat = (-90.0 + 0.5 * dy) + numpy.arange(nlat) * dy
    lon = 0.5 * dx + numpy.arange(nlon) * dx

    lon_coord = xarray.DataArray(
      dims=['lon'],
      coords=dict(lon=lon),
      data=lon,
      attrs=dict(standard_name='longitude', units='degrees_east')
    )

    lat_coord = xarray.DataArray(
      dims=['lat'],
      coords=dict(lat=lat),
      data=lat,
      attrs=dict(standard_name='latitude', units='degrees_north')
    )

    da = xarray.DataArray(
    name='land_area_fraction',
    data=land_frac,
    dims=["lat", "lon"],
    coords=dict(
        lat=lat_coord,
        lon=lon_coord,
    ),
    attrs=dict(
        STASH='m01s00i505',
        grid_staggering=6,
        valid_min=0.0,
        valid_max=1.0
        )
    )

    fracs_as_cube = da.to_iris()
    for dim in fracs_as_cube.coords():
        dim.bounds = None

    ants.io.save.netcdf(fracs_as_cube, out_fp + '/qrparm.landfrac')
    ants.io.save.ancil(fracs_as_cube, out_fp + '/qrparm.landfrac')

    # Also create the mask
    mask = da > 0.0
    mask.name = 'land_binary_mask'
    mask.attrs["STASH"] = 'm01s00i030'
    mask.attrs["grid_staggering"] = 6
    mask.attrs["valid_min"] = 0
    mask.attrs["valid_max"] = 1

    mask_as_cube = mask.to_iris()
    for dim in mask_as_cube.coords():
        dim.bounds = None

    ants.io.save.netcdf(mask_as_cube, out_fp + '/qrparm.mask')
    ants.io.save.ancil(mask_as_cube, out_fp + '/qrparm.mask')

if __name__ == '__main__':
    args = parse_args()
    ocean_mesh, ocean_mask = load_ocean_data(args.ocean_mesh_file)
    nlon, nlat = parse_resolution_string(args.atm_resolution)
    um_mesh = create_um_mesh(nlat, nlon)
    land_fractions = create_landfracs(um_mesh, ocean_mesh, ocean_mask)
    save_landfracs(nlat, nlon, land_fractions, args.output)
