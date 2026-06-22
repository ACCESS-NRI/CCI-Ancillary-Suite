# This simply sets all the land grid cells below -60 latitude to be 100% ice.
# The primary impetus for this change was the mismatch between what CCI
# considered land, and what the ocean model considered land. CCI treated the
# permanent ice shelves to be ocean e.g. Ross Ice Shelf, while the ACCESS ocean
# models treated it as land (in the form of ice).

# When ANTS attempts to fill in this mismatch using a spiral search, it was
# finding patches of bare ground and inland water on the CCI coast, and filling
# in the ice shelf with these. It's not clear why there was inland water at the
# antarctic coast, at locations where there should have been a permanent ice
# shelf.

import argparse
import ants
import numpy

def _parse_args():
    parser = ants.AntsArgParser(
            target_lsm=False,
            target_grid=False
            )

    parser.add_argument(
            '--ice-tile-id',
            type=int,
            default=9,
            help='Which tile ID is the model ice tile?'
            )

    return parser.parse_args()

def convert_to_ice(veg_frac, ice_tile_id):
    """
    Convert all the vegetation_area_fractions below -60 degrees latitude to be
    100% ice (where land).
    """

    # We want to create a mask which sets only the ice tile below -60 latitude to
    # 1.0, and every other land point below -60 latitude to 0.0.
    land_mask = ~veg_frac.data.mask
    lat_mask = veg_frac.coord('latitude').points > -60
    land_mask[:, lat_mask, :] = False
    
    ice_tile_mask = veg_frac.coord('pseudo_level').points.data == ice_tile_id

    ice_mask = ice_tile_mask[:, None, None] & land_mask
    nonice_mask = ~ice_tile_mask[:, None, None] & land_mask

    veg_frac.data[ice_mask] = 1.0
    veg_frac.data[nonice_mask] = 0.0
    return veg_frac

def main(veg_fractions_file, ice_tile, output, netcdf_only):
    vegetation_fractions = ants.load_cube(
            veg_fractions_file,
            'vegetation_area_fraction'
            )

    new_vegetation_fractions = convert_to_ice(vegetation_fractions, ice_tile)

    ants.io.save.netcdf(new_vegetation_fractions, output)
    if not netcdf_only:
        ants.io.save.ancil(new_vegetation_fractions, output)

if __name__ == '__main__':
    args = _parse_args()
    main(args.sources, args.ice_tile_id, args.output, args.netcdf_only)
