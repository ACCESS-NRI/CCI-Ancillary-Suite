import argparse
import numpy

import iris


def _parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
            '--tile-fractions',
            required=True,
            type=str,
            help='Tile fractions to apply minimum to.'
            )
    parser.add_argument(
            '--min-fraction',
            required=True,
            type=float,
            help='Minimum fraction to apply'
            )
    parser.add_argument(
            '--output',
            required=False,
            type=str,
            default='',
            help="""Output filename. If not supplied, then the --tile-fractions
            file is overwritten."""
            )
    parser.add_argument(
            '--netcdf-only',
            required=False,
            action='store_true',
            help="""Whether to only write output as NetCDF. If not specified,
            writes as both NetCDF and UM fields format."""
            )

    return parser.parse_args()


def apply_min_tile_fraction(area_fractions, min_fraction):
    """
    Remove all the tiles that have an area fraction less than min_fraction, and
    then re-normalise the area fractions so that the total fractions up to 0.0.
    """

    # The supplied land fractions should be a masked array, so we can apply a
    # simple set operation
    area_fractions.data[area_fractions.data < min_fraction] = 0.0

    # Now re-normalise
    new_totals = numpy.sum(area_fractions.data, axis=0)
    area_fractions.data = area_fractions.data / new_totals


if __name__ == '__main__':
    args = _parse_args()

    tile_fractions = iris.load_cube(args.tile_fractions,
                                    'vegetation_area_fraction')
    apply_min_tile_fraction(tile_fractions, args.min_fraction)


    # Note that it seems to be an ANTS "standard" that output file names are
    # supplied without a file extension, and then '.nc' is appended for NetCDF
    # files.

    # A string evaluates to false if it's length 0, so use that as the check
    # whether to write in place or to a new name
    if args.output:
        iris.fileformats.netcdf.save(tile_fractions, args.output + '.nc')

        if not args.netcdf_only:
            iris.fileformats.pp.save(tile_fractions, args.output)

    else:
        # Quickly do a naive check of what format the input file was- decide 
        # whether to drop ".nc" from the base name.
        if args.tile_fractions.endwith('.nc'):
            base_fname = args.tile_fractions[:-3]
        else:
            base_fname = args.tile_fractions

        iris.fileformats.netcdf.save(tile_fractions, base_fname + '.nc')

        if not args.netcdf_only:
            iris.fileformats.pp.save(tile_fractions, base_fname)
