import argparse
import ants
import numpy


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--land-fractions',
            required=True,
            help='Input land fractions used to create masks.'
            )

    parser.add_argument(
            '--output-path',
            required=True,
            help='Where to write generated land and sea masks.'
            )

    return parser.parse_args()


def landfrac_to_landmask(landfrac_path, out_path):
    land_fracs = ants.io.load.load_cube(landfrac_path)

    landmask = land_fracs.copy(data=(land_fracs.data > 0.0).astype(numpy.int64))
    landmask.attributes["valid_min"] = 0
    landmask.attributes["valid_max"] = 1
    landmask.attributes["grid_staggering"] = 6
    landmask.attributes["name"] = 'land mask'
    landmask.attributes["STASH"] = 'm01s00i30'

    for dim in landmask.coords():
        dim.bounds = None

    ants.io.save.netcdf(landmask, out_path + '/qrparm.mask')
    ants.io.save.ancil(landmask, out_path + '/qrparm.mask')

    seamask = land_fracs.copy(data=(land_fracs.data < 1.0).astype(numpy.int64))
    seamask.attributes["valid_min"] = 0
    seamask.attributes["valid_max"] = 1
    seamask.attributes["grid_staggering"] = 6
    seamask.attributes["name"] = 'sea mask'
    seamask.attributes["STASH"] = 'm01s00i30'

    for dim in seamask.coords():
        dim.bounds = None

    ants.io.save.netcdf(seamask, out_path + '/qrparm.mask_sea')
    ants.io.save.ancil(seamask, out_path + '/qrparm.mask_sea')


if __name__ == '__main__':
    args = _parse_args()
    
    landfrac_to_landmask(
            args.land_fractions,
            args.output_path,
            )
