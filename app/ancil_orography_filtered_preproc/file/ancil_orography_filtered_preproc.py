#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
"""
Filtered orography preprocessing
********************************

- Use an isotropic Raymond filter to filter features <= 6km of the provided
  surface altitude field.  See :func:`ants.analysis.filters.raymond`.
- Save to NetCDF (standard_name: 'surface_altitude').

"""
import ants
import ants.io.save as save
import iris
from ants.analysis.filters import raymond


def load_data(source):
    return ants.io.load.load_cube(source)


def main(sources, output, netcdf_only):
    source = load_data(sources)
    raymond(source, filter_length_scale=6000.0, isotropic=True)
    source.data[source.data < 0] = 0
    source.rename("surface_altitude_filtered")

    save.netcdf(source, output)
    if not netcdf_only:
        save.ancil(source, output)

    return source


if __name__ == "__main__":
    parser = ants.AntsArgParser()
    args = parser.parse_args()
    main(args.sources, args.output, args.netcdf_only)
