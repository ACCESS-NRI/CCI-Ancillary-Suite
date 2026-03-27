#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
"""
Unfiltered orography preprocessing
**********************************
- Save the source file to NetCDF (name: 'unfiltered_surface_altitude',
  STASH: 'm01s00i007'). This is to enable the data to be realised in
  sub-sets in the processor itself.

"""
import ants
import ants.fileformats
import ants.io.save as save
import iris


def load_data(source):
    return ants.io.load.load_cube(source)


def main(sources, output, netcdf_only):
    source = load_data(sources)
    source.rename("unfiltered_surface_altitude")
    source.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i007")

    save.netcdf(source, output)
    if not netcdf_only:
        save.ancil(source, output)

    return source


if __name__ == "__main__":
    parser = ants.AntsArgParser()
    args = parser.parse_args()
    main(args.sources, args.output, args.netcdf_only)
