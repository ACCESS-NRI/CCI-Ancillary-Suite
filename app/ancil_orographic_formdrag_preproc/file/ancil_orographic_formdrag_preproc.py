#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
"""
Orographic formdrag generation preprocessing
********************************************

The source file for orographic formdrag parameters has incorrect values set for 
grid staggering, lbcode and missing data. This preprocessor corrects those values.

In addition, the globe30 source data for the silhouette needs to be multipled by 
a factor of 1.823310 to account for the decrease in resolution in this source data. 
This is only needed for globe30 source data.

The following steps are taken:

* grid_staggering = 6
* bmdi = mule.REAL_MDI
* lbcode = 1
* multiply the field with STASH code of 17 by 1.823310

"""
import ants
import mule
from mule.operators import ScaleFactorOperator


def load_data(source):
    if len(source) != 1:
        raise RuntimeError(
            "Can only preprocess one file at a time, was given: {}".format(len(source))
        )
    source = source[0]
    ff = mule.AncilFile.from_file(source)
    return ff


def main(sources, output_filepath):
    scale_operator = ScaleFactorOperator(1.823310)
    ff = load_data(sources)
    ff.fixed_length_header.grid_staggering = 6
    for counter, field in enumerate(ff.fields):
        field.bmdi = mule._REAL_MDI
        field.lbcode = 1
        if field.lbuser4 == 17:
            ff.fields[counter] = scale_operator(field)
        elif field.lbuser4 == 18:
            continue
        else:
            raise ValueError(
                "Can only preprocess fields with STASH codes of 17 or 18, was given: {}".format(
                    field.lbuser4
                )
            )
    ff.to_file(output_filepath)
    return ff


if __name__ == "__main__":
    parser = ants.AntsArgParser()
    args = parser.parse_args()
    main(args.sources, args.output)
