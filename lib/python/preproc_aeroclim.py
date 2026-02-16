#!/usr/bin/env python
"""
Python script to pre-process netcdf GLOMAP MODE aerosol source files using
input PP files from a previous UM run.

"""

import argparse
import glob
import os
import warnings

import ants
import cf_units
import iris
import iris.fileformats.pp
import numpy as np

# Variable names used in LFRic interface to ancillary file
VAR_NAMES = {
    54103: "AiSo___103",
    54104: "AiSoSU_104",
    54105: "AiSoBC_105",
    54106: "AiSoOC_106",
    54107: "AcSo___107",
    54108: "AcSoSU_108",
    54109: "AcSoBC_109",
    54110: "AcSoOC_110",
    54111: "AcSoSS_111",
    54113: "CoSo___113",
    54114: "CoSoSU_114",
    54115: "CoSoBC_115",
    54116: "CoSoOC_116",
    54117: "CoSoSS_117",
    54119: "AiIn___119",
    54120: "AiInBC_120",
    54121: "AiInOC_121",
    54126: "NuSoOC_126",
}


def update_glomap_meta(cube, gregorian=False):
    """
    Routine to modify cube metadata to match what is required in GLOMAP MODE
    code in the UM and LFRic. This includes:
     * Adding a STASH code attribute that is the integer representation
     * Changing the variable name to match what is required in LFRic
     * Making sure the time coordinate matches what is required for the model
     * Added update_type and calendar_flexible attributes

    Parameters
    ----------
    cube : :class:`~iris.cube.Cube`
        Iris cube containing monthly mean data to be written out to file
    gregorian : boolean, optional
        Flag denoting calendar type. When true calendar is set to Gregorian, if
        false then set to 360day.

    """
    # Modify STASH section for field
    cube.attributes["STASH"] = iris.fileformats.pp.STASH(
        model=cube.attributes["STASH"].model,
        section=54,
        item=cube.attributes["STASH"].item,
    )
    # Required for GLOMAP MODE NetCDF files
    cube.attributes["stashcode"] = (
        cube.attributes["STASH"].section * 1000 + cube.attributes["STASH"].item
    )
    # Modify netcdf variable name
    cube.var_name = VAR_NAMES[cube.attributes["stashcode"]]
    # ANTS UKCA saver resets this as an integer but the UM expects a character
    #  string. Use of ncatted after the fact can convert it back to a string.
    cube.attributes["update_type"] = "2"
    # Time coordinate - convert to hours, taking account of any change in calendar
    date_bounds = cube.coord("time").units.num2date(cube.coord("time").bounds)
    calendar = cf_units.CALENDAR_STANDARD if gregorian else cf_units.CALENDAR_360_DAY
    time_units = cf_units.Unit("hours since 1970-01-01 00:00:00", calendar=calendar)
    cube.coord("time").units = time_units
    cube.coord("time").bounds = time_units.date2num(date_bounds)
    cube.coord("time").points = np.mean(cube.coord("time").bounds, axis=1)
    # This setting is not used in GLOMAP MODE climatologies at present
    cube.coord("time").attributes["calendar_flexible"] = "1"


def time_mean(cube, gregorian=False):
    """
    Routine to create mean on time coordinate using Cube.collapsed(). Also
    includes some metadata modification.

    Notes
    -----
    The time bounds for the complete climatology need to be strictly montonic
    and so the bounds are set to be the same as the first time point.
    Also, it is quicker for the data to be realised in this routine rather than
    leaving it as lazy and realising when the file is written out.

    Parameters
    ----------
    cube : :class:`~iris.cube.Cube`
        Iris cube to time mean
    gregorian : boolean, optional
        Flag denoting calendar type. When true calendar is set to Gregorian, if
        false then set to 360day.

    Returns
    -------
    : :class:`~iris.cube.Cube`
        Time mean with corrected metadata

    """
    # Create monthly mean
    mean = cube.collapsed("time", iris.analysis.MEAN)
    # Realise data - makes process a lot quicker!
    mean.data
    # Required for ANTS to accept file, needs to have monotonic time bounds?
    #  Also for UM/LFRic to be able to do date maths?
    mean.coord("time").points = cube[0].coord("time").points
    mean.coord("time").bounds = cube[0].coord("time").bounds
    # Update metadata
    update_glomap_meta(mean, gregorian=gregorian)
    return mean


def monthly_means(dpath, month, gregorian=False):
    """
    Wrapper routine to define and load input files, then produce the monthly
    mean climatology for all fields

    Notes
    -----
    Meaning all fields for a given month in turn is the most optimal route to
    generate a full climatology as file reads are much reduced

    Parameters
    ----------
    dpath : str
        Directory path containing input files
    month : str
        3 letter shorthand for month to mean over, defines files to load
    gregorian : boolean, optional
        Flag denoting calendar type. When true calendar is set to Gregorian, if
        false then set to 360day.

    Returns
    -------
    : :class:`~iris.cube.CubeList`
        Month climatology for each field and the given month

    """
    files = f"{dpath}/*{month}.pp"
    print(f"Reading {files}")
    # Load selected files using ANTS loader
    cubes = ants.io.load.load(sorted(glob.glob(files)))
    # Create monthly mean for each cube in cubes
    means = iris.cube.CubeList([time_mean(cube, gregorian=gregorian) for cube in cubes])
    return means


def parse_args():
    """
    Configure the argument parser for the script
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=str, required=True, help="Source root.")
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output filepath for pre-processed files",
    )
    parser.add_argument(
        "--calendar",
        type=str,
        choices=["gregorian", "360day"],
        default="gregorian",
        help="Calendar to use in output file, should be redundant if calendar-flexible were used in UM code",
    )
    parser.add_argument(
        "--datasource",
        type=str,
        required=False,
        help="Description of where source data has come from",
    )
    args = parser.parse_args()
    return args


def main():
    # Reduce output by ignoring expected and known warnings
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", message="Collapsing a non-contiguous coordinate")
    warnings.filterwarnings(
        "ignore", message="Unable to create instance of HybridHeightFactory"
    )

    # Read script arguments
    args = parse_args()

    # Turn calendar argument into a boolean flag
    gregorian = args.calendar == "gregorian"
    # Months to loop over
    months = [
        "jan",
        "feb",
        "mar",
        "apr",
        "may",
        "jun",
        "jul",
        "aug",
        "sep",
        "oct",
        "nov",
        "dec",
    ]
    # Loop over months to generate means for all fields in the given month
    clims = [
        monthly_means(args.source_root, month, gregorian=gregorian) for month in months
    ]
    # clims is a list of iris.cube.CubeList instances so unpack and merge to be
    #  a CubeList of monthly means for each field
    clims = iris.cube.CubeList([x for y in clims for x in y]).merge()
    # Add data source line to cubes for writing to NetCDF file.
    for clim in clims:
        clim.attributes["source"] = args.datasource
    print(clims)
    print(f"Writing {args.output}")
    # Save ancillary using UKCA saver
    ants.save(clims, args.output, saver="ukca")


if __name__ == "__main__":
    main()
