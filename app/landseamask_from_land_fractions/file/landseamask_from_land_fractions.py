import argparse
import iris


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--land-fractions',
            required=True,
            help='Input land fractions used to create masks.'
            )

    parser.add_argument(
            '--output',
            required=True,
            help='Where to write generated land and sea masks.'
            )

    return parser.parse_args()


def landfrac_to_landmask(landfrac_path, out_path):
    land_fracs = iris.load_cube(landfrac_path)

    landmask = land_fracs > 0.0
    seamask = land_fracs < 1.0

    landmask.attributes["valid_min"] = 0
    landmask.attributes["valid_max"] = 1
    landmask.attributes["grid_staggering"] = 6
    landmask.attributes["name"] = 'land mask'
    landmask.attributes["STASH"] = 'm01s00i30'

    for dim in landmask.coords():
        dim.bounds = None

    iris.fileformats.netcdf.save(landmask, out_path + 'qrparm.mask.nc')
    iris.fileformats.pp.save(landmask, out_path + 'qrparm.mask')

    seamask.attributes["valid_min"] = 0
    seamask.attributes["valid_max"] = 1
    seamask.attributes["grid_staggering"] = 6
    seamask.attributes["name"] = 'sea mask'
    seamask.attributes["STASH"] = 'm01s00i30'

    for dim in seamask.coords():
        dim.bounds = None

    iris.fileformats.netcdf.save(seamask, out_path + 'qrparm.mask_sea.nc')
    iris.fileformats.pp.save(seamask, out_path + 'qrparm.mask_sea')



if __name__ == '__main__':
    args = _parse_args()
    
    landfrac_to_masks(
            args.land_fractions,
            args.output_landmask,
            args.output_seamask
            )
