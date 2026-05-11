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

def _parse_args():
    parser = argparse.ArgumentParser("""Script to set all vegetation below -60
            latitude to permanent ice.""")

    parser.add_argument(
            '-i',
            '--input',
            type=str,
            required=True,
            help='Target vegetation_area_fractions to convert.'
            )

    parser.add_argument(
            '--ice-tile-id',
            type=int,
            default=9,
            help='Which tile ID is the model ice tile?'
            )

    parser.add_argument(
            '-o',
            '--output',
            type=str,
            required=True,
            help='Filename to write to.'
            )

    return parser.parse_args()

def convert_to_ice(veg_frac, ice_tile_id):
    """
    Convert all the vegetation_area_fractions below -60 degrees latitude to be
    100% ice (where land).
    """

    # ANTS uses a float NaN of very large- so filter by only selecting lower
    # than a value. Take the original mask, which masks out the ocean, and then
    # add the latitude and tile specific mask to allow masked indexing.
    ice_mask = c.data.mask.copy
    lat_mask = veg_frac.coord('latitude').points > -60
    ice_tile = veg_frac.coord('pseudo_level').points == ice_tile_id
    ice_mask[~ice_tile, lat_mask, :] = True
    nonice_mask[ice_tile, lat_mask, :] = True
    
    veg_frac.data[~ice_mask] = 1.0
    veg_frac.data[~nonice_mask] = 0.0
    return veg_frac

def main(veg_fractions_file, ice_tile, output):
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
    main(args.input, args.ice_tile, args.output)
