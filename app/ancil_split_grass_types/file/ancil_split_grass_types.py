import ants
import json
import numpy
import scipy

# This script uses an external dataset to split a generic grass type from a
# source land fractions using an external dataset. The process is split into 3
# stages:
# 1. Process a provided mapping to determine which tiles in the external lan
#   cover dataset map to which tiles in the source land fractions. The mapping
#   is provided as a dictionary of {source_tiles: [external_tiles]}, e.g.:
#   
#   {6: [7, 8], 7: [9, 10], 8: [11]}
#   
#   Would be a mapping that uses tiles 7, 8 from the external dataset to inform
#   tile 6 from the source land cover, 9, 10 to inform tile 7, 11 to inform 8.
# 2. Use the mapping to reduce the external dataset to an array of size
#   (len(mapping), nlat, nlon), where the first index describes the fraction of
#   grass that should be assigned to each of the grass tiles in the source
#   land fractions.
# 3. Use the array to override the original grass fractions in the source land
#   fractions dataset.

def _get_parser():
    parser = ants.AntsArgParser()
    parser.add_argument(
            '--original-fractions',
            type=str,
            required=True,
            help='Original land fractions to split grasses on.'
            )
    parser.add_argument(
            '--grass-mapping',
            type=str,
            required=True,
            help='JSON file describing the source and target grass tiles.'
            )

    return parser

def process_grass_mapping(grass_mapping_file):
    """
    Read the grass mapping JSON and check the contents are valid. Casts all the
    keys to integers, and values to lists of integers that can be used for
    indexing.
    """

    with open(grass_mapping_file, 'r') as f:
        grass_mapping = json.load(f)

    # We want to rewrite the keys and values to be integers- make a copy dict
    int_grass_mapping = {}
    for k, v in grass_mapping.items():
        try:
            k = int(k) - 1
        except:
            ValueError(f'The mapping key {k} could not be cast to an int')

        if not isinstance(v, list):
            v = [v]

        for i in range(len(v)):
            try:
                v[i] = int(v[i]) - 1
            except:
                ValueError(f'The value {v} could not be cast to an int')

        int_grass_mapping[k] = v

    return int_grass_mapping
        
def load_data(grass_source, orig_fractions):
    """
    Load the datasets for the grass fractions and the original land cover.
    """
    grasses = ants.io.load.load_cube(grass_source)
    original_fractions = ants.io.load.load_cube(orig_fractions)

    return grasses, original_fractions

def preprocess_grass_source(grasses, grass_mapping):
    """
    Take the indices specified by the values in the grass mapping, and
    consolidate them into single pages of an array for each key in the mapping.
    """

    # Prepare the condensed fractions- same land shape, with one page for each
    # entry in the mapping dict.
    condensed_data = numpy.ndarray(
        (len(grass_mapping), *grasses.data.shape[1:]),
        dtype=float
        )

    for target_ind, (_, source_inds) in enumerate(grass_mapping.items()):
        condensed_data[target_ind, :, :] = numpy.sum(grasses.data[source_inds, :, :], axis=0)

    # Now normalise to that the fractions along the tile axis sum to 1.0
    net_fractions = numpy.sum(condensed_data, axis=0)

    condensed_data /= net_fractions

    # This will result in a nan wherever the NCAR thought there was no grass.
    # We don't want to apply this to the original fractions- the assumption we
    # made was that CCI will determine how much grass there is, and NCAR splits
    # it, so we should keep to that rule. So we will fill in the NaN entries
    # that resulted from zero fractions using nearest neighbour.
    mask = numpy.where(~numpy.isnan(condensed_data))
    interp = scipy.interpolate.NearestNDInterpolator(mask, condensed_data[mask])
    condensed_data = interp(*numpy.indices(condensed_data.shape))

    return condensed_data

def apply_grass_fractions(orig_fractions, grass_fractions, grass_mapping):
    """
    Apply the processed grass fractions to the grass types in the original
    land fractions.
    """

    # First, determine the total grass fractions in each cell in the original
    # fractions
    orig_grass_ids = list(grass_mapping.keys())
    orig_grasses = numpy.sum(orig_fractions.data[orig_grass_ids, :, :], axis=0)

    for source_ind, target_ind in enumerate(orig_grass_ids):
        orig_fractions.data[target_ind, :, :] = grass_fractions[source_ind, :, :] * orig_grasses

    # We also need to think about what happens when the CCI fractions say there
    # is grass, but the NCAR fractions say there is none. At this stage, it
    # will end up removing all the grass without compensating the other types.
    return orig_fractions

def main(grass_source, original_fractions, mapping_file, output, netcdf_only):
    """
    Reprocess the grass fractions in the original fractions.
    """
    grass_mapping = process_grass_mapping(mapping_file)
    print(f'The generated grass mapping: {grass_mapping}')
    grass_source, orig_fractions = load_data(grass_source, original_fractions)
    
    grass_fractions = preprocess_grass_source(grass_source, grass_mapping)
    print(f'The shape of the processed grass fractions: {grass_fractions.shape}')
    new_fractions = apply_grass_fractions(
        orig_fractions,
        grass_fractions,
        grass_mapping
        )

    ants.io.save.netcdf(new_fractions, output)
    if not netcdf_only:
        ants.io.save.ancil(new_fractions, output)

if __name__ == '__main__':
    args = _get_parser().parse_args()
    main(
            args.sources,
            args.original_fractions,
            args.grass_mapping,
            args.output,
            args.netcdf_only
            )
