import xarray

def parse_args():
    parser = argparse.ArgumentParser(
            """Script to merge the CCI land cover for a specific year with the
            permanent water bodies dataset. The land cover dataset is at 300m
            resolution, while the permanent water is at 150m resolution."""
            )
    parser.add_argument(
            '--land-cover',
            type=str,
            help='Source land cover dataset',
            required=True
            )
    parser.add_argument(
            '--water-bodies',
            type=str,
            help='Source water bodies dataset',
            required=True
            )
    parser.add_argument(
            '-o',
            '--output',
            help='Output path',
            required=True
            )

    return parser.parse_args()


def main(land_cover_source, water_bodies_source, output):
    # Load in the data sources
    land_cover = xarray.open_dataset(land_cover_source)
    water_bodies = xarray.open_dataset(water_bodies_source)

    upscaled_water_bodies = upscale_water_bodies(water_bodies)

    overlay_water_on_land_cover(upscaled_water_bodies, land_cover)

    land_cover.to_netcdf(output)

def upscale_water_bodies(water_bodies):
    """
    Upscale the provided water bodies dataset at 150m resolution to 300m
    resolution. The water bodies dataset is a categorical dataset, with
    0=ocean, 1=land, 2=inland_water.
    """

    # We should be able to do this somewhat cleverly, instead of relying on
    # some windowed most common algorithm. Given we know the only possible
    # values are 0, 1, 2, we apply a three step process:
    # 1. Change the value "2" to "4" everywhere
    # 2. Perform the windowed sum, so we get the sum of values in every 2x2 box
    # 3. Inspect the resulting value, to determine what the most common value
    #   was in the original 2x2 box.
    # The bias will point "downwards"; the lowest flag value will take
    # precedence when two flags are equally most common.
    # This means the rules when inspecting the resulting values are:
    # 1. If value // 4 > 2, then the most common flag was a 4 (originally 2)
    # 2. If mod(value, 4) > 2, then the most common flag was a 1
    # 3. Otherwise, the most common flag was a 0.
    water_data = water_bodies["water_classification"].data
    water_data[water_data == 2] = 4

    upscaled_data = water_data[1:2:, 1:2:] + water_data[2:2:, 1:2:] + \
                    water_data[1:2:, 2:2:] + water_data[2:2:, 2:2:]

    # Initialise new arrays of 0s of same shape, so we can just override the
    # cells that need changing.
    flag_data = numpy.zeros_like(upscaled_data)
    flag_data[upscaled_data // 4 > 2] = 2
    flag_data[numpy.mod(upscaled_data, 4) > 2] = 1
    
    return flag_data

def overlay_water_on_land_cover(water, land_cover):
    """
    Apply the permanent water map over the top of the land cover map. The
    permanent water map is taken as the absolute truth in regards to water, and
    will override any conflicting land values in the land cover map. Any cells
    considered water by the land cover map, but land by the permanent water
    map, will be filled by the nearest valid land point i.e. nearest point
    along the same latitude considered land by the permanent water map.
    """
    
    flag_meanings = land_cover['lccs_class'].attrs['flag_meanings'].split()
    flag_values = land_cover['lccs_class'].attrs['flag_values']

    # Identify the flag associated with water
    water_index = flag_meanings.index('water')
    water_flag = flag_values[water_index]

    # Split the water into sea_ocean_water and inland_water. To make it clear
    # that this is a different flag to the original water, don't reuse the same
    # flag- use the next two free flag values
    sea_ocean_water_flag = water_flag
    while sea_ocean_water_flag in land_cover['lccs_class'].attrs['flag_values']:
        sea_ocean_water_flag += 1

    inland_water_flag = water_flag
    while inland_water_flag in land_cover['lccs_class'].attrs['flag_values']:
        inland_water_flag += 1

    # Now insert them in the locations adjacent to the original water flag
    flag_meanings.insert(water_index+1, 'sea_ocean_water')
    flag_meanings.insert(water_index+2, 'inland_water')
    flag_values.insert(water_index+1, sea_ocean_water_flag)
    flag_values.insert(water_index+2, inland_water_flag)

    # Remove the now-unused water flag
    flag_meanings.pop(water_index)
    flag_values.pop(water_index)

    # First, set all original water values to N/A
    new_land_cover = xarray.where(land_cover['lccs_class'] != water_flag)
    new_land_cover[water['water_classification'] == 0] == sea_ocean_water_flag
    new_land_Cover[water['water_classification'] == 2] == inland_water_flag

    # Now the only remaining N/A values are the ones that were water according
    # to the original land cover, but considered land by the permanent water.
    new_land_cover.interpolate_na(dim='lon', method='nearest')
    land_cover['lccs_class'] = new_land_cover

def load_data(sources):
    # Expect sources to be a list of two strings, describing the land cover
    # and water body datasets
    land_cover = ants.io.load.load_cube(sources[0], "land_cover_lccs")
    water_bodies = ants.io.load.load_cube(sources[1], "water_classification")

    return land_cover, water_bodies

def 
if __name__ == "__main__":
    args = _parse_args()

    # The sources contains [CCI_land_cover, CCI_water_bodies]
    main(args.land_cover, args.water_bodies, args.output)

