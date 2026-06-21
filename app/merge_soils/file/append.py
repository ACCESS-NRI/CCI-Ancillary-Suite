#!/usr/bin/env python
"""
The CAP reads specifically 10 fields from the input file when generating the
snow ancillaries (ancilSmcsnow).  There is more than 10 fields present in the
soils ancillary which means that ancilSmcsnow is dependent on field order...
We move those fields it doesn't need to the end of the file to avoid a
segmentation fault.

"""
import warnings

import ants
import ants.io.save as save
import iris

iris.FUTURE.save_split_attrs = True

def main(filenames, output):
    cubes = ants.io.load.load(filenames)

    # Remove scalar time coordinate if it exists in any of the cubes
    # to ensure that all cubes have the smae coordinates
    for i,cube in enumerate(cubes):
        if "time" in [coord.name() for coord in cubes[i].aux_coords]:
            cubes[i].remove_coord("time")

    # Remove history attribute
    for cube in cubes:
        cube.attributes.pop('history')

    # Sort by the following stash items (lbuser4).
    items = [40, 41, 43, 207, 47, 44, 46, 48, 220, 223, 8, 418, 419, 420]
    end = len(items) + 1
    
    def sortme(cube):
        stash = cube.attributes['STASH']
        try:
            order = items.index(stash.item)
        except ValueError:
            warnings.warn('Was this stash expected?? {}'.format(stash))
        return order
    cubes = sorted(cubes, key=sortme)
    save.ancil(cubes, output)
    save.netcdf(cubes, output)


if __name__ == '__main__':
    parser = ants.AntsArgParser()
    args = parser.parse_args()
    main(args.sources, args.output)

