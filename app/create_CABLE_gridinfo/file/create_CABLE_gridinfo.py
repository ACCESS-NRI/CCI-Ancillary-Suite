import xarray

# The mapping from CABLE gridinfo variables to JULES ancillaries:
# isoil => no equivalent variable
# SoilMoist => state variable, not ancillary.
# SoilTemp => state_variable, not ancillary.
# SnowDepth => state_variable, not ancillary.
# Albedo => Unused, set to 0.2 across bands
# LAI => qrparm.veg.func["leaf_area_index"]
# SoilOrder => no equivalent variable
# Ndep => unclear
# Nfix => unclear
# Pwea => unclear
# Pdust => unclear
# clay => qrparm.soil.nc["T_CLAY"], no unit conversion
# silt => qrparm.soil.nc["T_SILT"], no unit conversion
# sand => qrparm.soil.nc["T_SAND"], no unit conversion
# swilt => qrparm.soil.nc["sm_wilt"], 
# sfc => qrparm.soil.nc["sm_crit"]
# ssat => qrparm.soil.nc["sm_sat"]
# bch => unclear
# hyds => qrparm.soil.nc["satcon"]
# sucs => qrparm.soil.nc["sathh"]
# rhosoil => qrparm.soil.nc["BD"]
# cnsd => qrparm.soil.nc["hcon"]
# css => qrparm.soil.nc["hcap"]
# albedo2 => qrparm.soil.nc["soil_albedo"]
def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            "--veg-frac",
            type=str,
            required=True,
            description="File containing the vegetation fractions."
            )
    
    parser.add_argument(
            "--veg-func",
            type=str,
            required=True,
            description="File containing the LAI and canopy heights."
            )

    parser.add_argument(
            "--soil-properties",
            type=str,
            required=True,
            description="File containing soil properties."
            )

