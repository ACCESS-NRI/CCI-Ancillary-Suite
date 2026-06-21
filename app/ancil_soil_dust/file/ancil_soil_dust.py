#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
"""
Soil dust fields
****************

Soil dust fields represents the relative mass fraction of soil particles in
each of the 6 size bins - the term `M_rel` in this 2001 paper (https://doi.org/10.1029/2000JD900795).

The obs-based soil dataset contains only the clay, silt and sand fractions (by
mass) of the soil, i.e. the fraction of the total soil mass which has particles
within the size ranges for "clay", "silt" and "sand". In order to relate these
values to the fraction of soil mass in each of the bins in the dust parent soil
ancillary we have to make assumptions about the size distributions of each.
This is done in `dM/dlog(R)` vs `log(R)` space (as in fig.1 of the 2001 paper
- the x-axis there is labelled `(R)`, but you'll see that it is on a log
scale).

Note that `log(R)` is used to mean `log10(R)` here.  The clay, silt and sand
fractions in the soil dataset essentially give the integral of the `dM/dlog(R)`
distribution across the appropriate R ranges for each - the "areas under the
curves" in fig.1  By defining the form of the soil size distribution, we can
use the measured data to obtain the actual size distribution and then by
integrating across the limits of the model bins, we can calculate the fraction
of the mass within each.

Fig.1 illustrates what assumptions are made about the size distributions of the
clay, silt and sand fractions: the boundary between clay and silt is taken to
be at 1 micro meter radius, and that between silt and sand at 25 micro meters;
the maximum radius of sand is taken to be 2000 micro meters (the 200 micro
meters limit in fig.1 is a mistake that slipped through). For silt and sand, a
simple bin distribution (`dM/dlog(R)` constant across the range of `log(R)`) is
assumed for each, giving the "histogram-like" shapes for silt and sand shown.
The situation is a bit more complex for clay as there is no minimum size of
measured clay particles. Here we assume that `dM/dlog(R)` is a linear function
of `log(R)`, with continuity at the clay-silt boundary
(`dm/dlog(R) = C log(R) + D`, where `C` and `D` are constant; and at
`R=1.0 micro meters` `dM_clay/dlog(R)_clay = dM_silt/dlog(R)`). For the given
clay and silt fractions, this defines the slope of the `dM/dlog(R)` line for
clay and the minimum clay size (`A` in fig.1). From here, the clay, silt and
sand fractions allow us to calculate the `dM/dlog(R)` distribution for all
particle sizes. Then we can integrate this across each of the bins to obtain
the mass fraction in each bin.

Contact:
Stephanie Woodward (stephanie.woodward@metoffice.gov.uk)

"""
import ants
import ants.io.save as save
import iris
import numba
import numpy as np

iris.FUTURE.save_split_attrs = True

# max-min particle dims for each division.
MAX_DIM_PART = np.array([.2e-6, .632456e-6, 2e-6, 6.32456e-6, 2e-5,
                         6.32456e-5])
MIN_DIM_PART = np.array([.0632456e-6, .2e-6, .632456e-6, 2e-6, 6.32456e-6,
                         2e-5])
# max sand-silt-clay particle diameter.
MAX_DIM_CLAY = 2e-6
MAX_DIM_SILT = 50e-6
MAX_DIM_SAND = 2000e-6

# minimum silt/clay fraction for dust calculation.
MIN_FRACTION = 1e-8

# Number of soil fraction divisions
NSFPDS = 10

# Cube metadata for our soil particles for each division.
METADATA = [dict(long_name=f'mass_fraction_of_soil_particles_in_ukmo_division{i}', units='kg/kg', STASH=f'm01s00i42{i}') for i in range(1, 6+1)]


@numba.jit(nopython=True, nogil=True)
def _calc_soil_dust_mass_fraction(f_sand, f_silt, f_clay, land_sea_targ,
                                  landice_mask, max_dim_part, min_dim_part):
    """
    Calculate the soil dust mass fraction for each soil particle division.

    Calculated on a point by point basis so works with 1D arrays.

    Parameters
    ----------
    f_sand : np.ndarray
        1D fraction of sand.
    f_silt : np.ndarray
        1D fraction of silt.
    f_clay : np.ndarray
        1D fraction of clay.
    landice_mask : np.ndarray
        1D Land ice mask, where True denotes ice while False denotes no ice.
    land_sea_tar : np.ndarray
        1D Land sea mask, where True denotes land and False denotes ocean.
    max_dim_part : np.ndarray
        Max particle dimensions for each division.
    min_dim_part : np.ndarray
        Min particle dimensions for each division.

    Returns
    -------
    : np.ndarray
        Returns a numpy array of shape [ndiv, npts] representing the soil dust
        mass fraction of soil particles in each ukmo division.

    """
    # Number of soil divisions
    numsubs = max_dim_part.size
    p = np.zeros((numsubs,) + f_silt.shape)

    beta = np.log10(MAX_DIM_CLAY)
    for i in range(f_clay.shape[0]):
        if not land_sea_targ[i] or landice_mask[i]:
            # Skip over ocean and ice points.
            continue

        # Apply minimum silt and clay fraction to aid soil dust calculation.
        if f_silt[i] < MIN_FRACTION:
            f_silt[i] = MIN_FRACTION
        if f_clay[i] < MIN_FRACTION:
            f_clay[i] = MIN_FRACTION

        csi = f_silt[i] / (np.log10(MAX_DIM_SILT) -
                           np.log10(MAX_DIM_CLAY))
        alpha = beta - (2.0*(f_clay[i]/csi))

        for numsub in range(numsubs):
            log_size_sub = (np.log10(max_dim_part[numsub]) -
                            np.log10(min_dim_part[numsub])) / NSFPDS
            size_sub = 10**log_size_sub
            mass_division = 0.

            for nsfpd in range(NSFPDS):
                log_min_dim_sub = (np.log10(min_dim_part[numsub]) +
                                   (nsfpd * log_size_sub))
                min_dim_sub = 10**log_min_dim_sub
                log_max_dim_sub = log_min_dim_sub + log_size_sub
                max_dim_sub = 10**log_max_dim_sub
                log_rep_size = (np.log10(max_dim_sub) +
                                np.log10(min_dim_sub))/2
                rep_size = 10**((np.log10(max_dim_sub) +
                                np.log10(min_dim_sub))/2)

                if log_max_dim_sub <= alpha:
                    pp = 0

                elif log_min_dim_sub <= alpha:
                    pp = (csi*(log_max_dim_sub - alpha) *
                          ((csi*(log_max_dim_sub + alpha))/2. +
                           (2.*f_clay[i]) - (beta*csi)) /
                          (2.*f_clay[i]))
                elif max_dim_sub <= MAX_DIM_CLAY:
                    pp = csi*log_size_sub*((csi*log_rep_size +
                                            2*f_clay[i] - beta*csi) /
                                           (2*f_clay[i]))

                elif min_dim_sub <= MAX_DIM_CLAY:
                    pp = (((csi*(np.log10(MAX_DIM_CLAY) -
                                 log_min_dim_sub)) /
                           (2*f_clay[i])) *
                          (((csi*(np.log10(MAX_DIM_CLAY) +
                                  log_min_dim_sub)) /
                            2) +
                           (2*f_clay[i]) -
                           (np.log10(MAX_DIM_CLAY) * csi)) +
                          csi*(log_max_dim_sub - np.log10(MAX_DIM_CLAY)))

                elif max_dim_sub <= MAX_DIM_SILT:
                    pp = csi*log_size_sub

                elif min_dim_sub <= MAX_DIM_SILT:
                    pp = (csi*(np.log10(MAX_DIM_SILT) -
                               np.log10(min_dim_sub)) +
                          f_sand[i] * (np.log10(max_dim_sub) -
                                       np.log10(MAX_DIM_SILT)) /
                          (np.log10(MAX_DIM_SAND) - np.log10(MAX_DIM_SILT))
                          )

                elif max_dim_sub <= (MAX_DIM_SAND * (1.+1.e-5)):
                    pp = ((f_sand[i] * log_size_sub) /
                          (np.log10(MAX_DIM_SAND) - np.log10(MAX_DIM_SILT))
                          )

                mass_division += abs(pp)
            p[numsub, i] = mass_division
    return p


def calc_soil_dust_mass_fraction(f_sand, f_silt, f_clay, land_sea_targ,
                                 landice_mask):
    """
    Calculate soil dust mass fraction of soil particles in each ukmo division.

    Parameters
    ----------
    f_sand : np.ndarray
        Fraction of sand.
    f_silt : np.ndarray
        Fraction of silt.
    f_clay : np.ndarray
        Fraction of clay.
    landice_mask : np.ndarray
        Land ice mask, where True denotes ice while False denotes no ice.
    land_sea_tar : np.ndarray
        Land sea mask, where True denotes land and False denotes ocean.

    Returns
    -------
    : np.ndarray
        Mass fraction of soil particles in each ukmo division.

    """
    # Check that divisions cover whole range up to max value for sand.  If not,
    # add an extra division which isn't written out.
    max_dim_part = MAX_DIM_PART
    min_dim_part = MIN_DIM_PART
    if MAX_DIM_PART[-1] < MAX_DIM_SAND:
        max_dim_part = np.hstack([MAX_DIM_PART, [MAX_DIM_SAND]])
        min_dim_part = np.hstack([MIN_DIM_PART, [MIN_DIM_PART[-1]]])
    return _calc_soil_dust_mass_fraction(f_sand.reshape(-1), f_silt.reshape(-1),
                                         f_clay.reshape(-1), land_sea_targ.reshape(-1),
                                         landice_mask.reshape(-1), max_dim_part,
                                         min_dim_part)


def calc_soil_dust_fields(sources):
    """
    Calculate soil dust mass fraction of soil particles in each ukmo division.

    Parameters
    ----------
    sources : :class:`iris.cube.CubeList`
        CubeList containing the sand, silt and clay fractions as well as the
        land cover type fraction field (from which we extract the landsea mask
        and ice fraction fields).

    Returns
    -------
    : :class:`iris.cube.CubeList`
        Mass fraction of soil particles in each ukmo division.

    """
    sand, silt, clay, lctf = sources

    landice_mask = lctf[-1].data
    landsea_mask = ~np.ma.getmaskarray(landice_mask)
    # NOTE: We might want in future to check that lctf and sand, silt, clay
    # masks are all the same.

    sand_data = np.asarray(sand.data)
    silt_data = np.asarray(silt.data)
    clay_data = np.asarray(clay.data)
    landice_mask = np.asarray(landice_mask.data)
    landice_mask = (landice_mask > 0) * landsea_mask

    res = calc_soil_dust_mass_fraction(sand_data, silt_data, clay_data, landsea_mask,
                                       landice_mask)

    # Turn resulting arrays into cubes.
    cubes = []
    for i in range(6):
        cube = sand.copy()
        cube.rename(METADATA[i]['long_name'])
        cube.units = METADATA[i]['units']
        cube.attributes['STASH'] = iris.fileformats.pp.STASH.from_msi(METADATA[i]['STASH'])
        cube.data = np.ma.array(res[i].reshape(cube.shape), mask=~landsea_mask)
        cubes.append(cube)

    return cubes


def load(sand_silt_clay_frac_path, lctf_path):
    sand_con = iris.AttributeConstraint(STASH='m01s00i420')
    silt_con = iris.AttributeConstraint(STASH='m01s00i419')
    clay_con = iris.AttributeConstraint(STASH='m01s00i418')
    sand, silt, clay = ants.io.load.load_cubes(sand_silt_clay_frac_path,
                                       [sand_con, silt_con, clay_con])
    lctf_stash_con = iris.AttributeConstraint(STASH='m01s00i216')
    lctf = ants.io.load.load_cube(lctf_path, lctf_stash_con)
    return sand, silt, clay, lctf


def main(sand_silt_clay_frac_path, lctf_path, output, netcdf_only):
    sand, silt, clay, lctf = load(sand_silt_clay_frac_path, lctf_path)
    #res = decomp.decompose(calc_soil_dust_fields, [sand, silt, clay, lctf])
    res = calc_soil_dust_fields([sand, silt, clay, lctf])

    # Append the sand, silt and clay to the soil dust file we generate.
    res.extend([sand, silt, clay])

    save.netcdf(res, output)
    if not netcdf_only:
        save.ancil(res, output)


if __name__ == '__main__':
    parser = ants.AntsArgParser()
    lct_help = (
        "Land cover type fraction ancillary, used for ensuring "
        "consistency of the topographic index with the ice field and "
        "the landsea mask."
    )
    parser.add_argument("--lct-ancillary", type=str, help=lct_help, required=True)
    args = parser.parse_args()
    main(args.sources, args.lct_ancillary, args.output, args.netcdf_only)
