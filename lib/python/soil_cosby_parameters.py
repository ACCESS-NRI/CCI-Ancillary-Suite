#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
'''
Calculates soil parameters from the soil textural class or from the fraction
of sand, silt and clay in the mineral soil component.
The parameters can also be adjusted for a organic soil component.

 - Karina Williams and Heather Ashton

'''

import argparse
import numpy as np

'''
References
  * [Cosby 1984]
    Cosby, B.J., Hornberger, G.M., Clapp, R.B. and Ginn, T.T., 1984. A
    statistical exploration of the relationships of soil moisture
    characteristics to the physical properties of soils. Water Resources
    Research. 20. 682-690.
    http://onlinelibrary.wiley.com/doi/10.1029/WR020i006p00682/pdf
    Note that logs in this paper are log to the base 10.
  * [AncilDoc_Cosby]
    AncilDoc_Cosby.html in the ancil documentation (see attachment to
    https://code.metoffice.gov.uk/trac/ancil/ticket/250)
  * [AncilDoc_SoilThermal]
    AncilDoc_SoilThermal.html in the ancil documentation (see attachment to
    https://code.metoffice.gov.uk/trac/ancil/ticket/250)
  * [Dankers 2011]
    R. Dankers, E. J. Burke, J. Price, Simulation of permafrost and seasonal
    thaw depth in the JULES land surface scheme
    www.the-cryosphere.net/5/773/2011/
  * [Chadburn 2015]
    S. Chadburn, E. Burke, R. Essery, J. Boike, M. Langer, M. Heikenfeld,
    P. Cox, and P. Friedlingstein, An improved representation of physical
    permafrost dynamics in the JULES land surface model
    www.geosci-model-dev-discuss.net/8/715/2015/
  * [Best 2011]
    M. J. Best, M. Pryor, D. B. Clark, G. G. Rooney, Essery, C. B. Menard, 
    J. M. Edwards, M. A. Hendry, A. Porson, N. Gedney, L. M. Mercado, S. Sitch, 
    E. Blyth, O. Boucher, P. M. Cox, C. S. B. Grimmond, R. J. Harding
    The Joint UK Land Environment Simulator (JULES), model description -
    Part 1: Energy and water fluxes
    http://www.geosci-model-dev.net/4/677/2011/
'''


# array contains percentage silt, percentage sand, percentage clay from
# table 2 of [Cosby 1984]
# plus silt (see email from Heather Ashton 6.7.15)
SOIL_TEXTURAL_CLASS = {
    'Sand':            [ 5, 92,  3],
    'Loamy sand':      [12, 82,  6],
    'Sandy loam':      [32, 58, 10],
    'Loam':            [39, 43, 18],
    'Silty loam':      [70, 17, 13],
    'Sandy clay loam': [15, 58, 27],
    'Clay loam':       [34, 32, 34],
    'Silty clay loam': [56, 10, 34],
    'Sandy clay':      [ 6, 52, 42],
    'Silty clay':      [47,  6, 47],
    'Clay':            [20, 22, 58],
    'Silt':            [90,  7,  3]
}


# Parameters for organic soil used in the SOC experiment
# from [Dankers 2011] table 2. Layers are 0-10cm, 10-35cm, 35-100cm
ORGANIC_SOIL = {'b': [2.7, 6.1, 12.0],
                # b, Brooks-Corey exponent in the soil hydraulic
                # characteristics
                'sathh': [0.0103, 0.0102, 0.0101],
                # Psi_s, absolute value of the soil matric suction at
                # saturation, m
                'satcon': [0.28, 0.002, 0.0001],
                # K_s, hydraulic conductivity at saturation, kg m-2 s-1
                'sm_sat': [0.93, 0.88, 0.83],
                # theta_s, volumetric soil moisture content at saturation,
                # m3 water per m3 soil
                'sm_crit': [0.11, 0.34, 0.51],
                # theta_c, volumetric soil moisture content at the critical
                # point at which soil moisture stress starts to reduce
                # transpiration, m3 water per m3 soil
                'sm_wilt': [0.03, 0.18, 0.37],
                # theta_w, volumetric soil moisture content at the wilting
                # point, m3 water per m3 soil
                'hcap': [0.58E+6, 0.58E+6, 0.58E+6],
                # c, dry heat capacity, J m-3 K-1
                'hcon': [0.06, 0.06, 0.06]
                # lambda, dry thermal conductivity, W m-1 K-1
                }


COSBY_PARAMETER_DESCRIPTION = {
                'b': 'Brooks-Corey exponent in the soil hydraulic characteristics, dimensionless',
                'sathh': 'absolute value of the soil matric suction at saturation (Psi_s), in m',
                'satcon': 'hydraulic conductivity at saturation (K_s), in kg m-2 s-1',
                'sm_sat': 'volumetric soil moisture content at saturation (theta_s), in m3 water per m3 soil',
                'sm_crit': 'volumetric soil moisture content at the critical point at '
                           'which soil moisture stress starts to reduce transpiration (sm_crit), in m3 water per m3 soil',
                'sm_wilt': 'volumetric soil moisture content at the wilting point (sm_wilt), in m3 water per m3 soil',
                'hcap': 'dry heat capacity (c), in J m-3 K-1',
                'hcon': 'dry thermal conductivity (lambda), in W m-1 K-1',
                'albsoil': 'soil albedo (not calculated by this script)',
                'BD': 'bulk density of the soil in kg per m3',
                'soil_carb': 'amount of soil carbon in the soil in kg of carbon per m3 of soil',
                'clay': 'soil clay fraction (g clay per g soil)',
                }


RHO_W = 999.97  # density of water kg m-3 (at 4 degC)
G_GRAVITY = 9.81  # acceleration due to gravity in m s-2


class _SoilType(object):
    PARAMETERS = ['b', 'sathh', 'satcon', 'sm_sat', 'sm_crit', 'sm_wilt',
                  'hcap', 'hcon', 'albsoil', 'clay']

    ''' defines the properties of the soil '''
    def __init__(self, f_clay, f_sand, f_silt, bulk_den, org_carb,
                 f_organic=None, soc_layer=None):
        '''
        Kwargs:

         - f_clay : fraction of clay in mineral soil
         - f_sand : fraction of sand in mineral soil
         - f_silt : fraction of silt in mineral soil
         - bulk_den  : Bulk Density (kg / nm3)
         - org_carb  : Organic carbon (% weight)
         - f_organic : fraction of organic matter in soil
         - soc_layer : layer to take organic properties from (described in
           [Dankers 2011] table 2)

        '''
        total_fraction = sum([f_clay, f_sand, f_silt])
        if np.max(abs(1.0 - total_fraction)) > 1.0E-3:
            raise ValueError('Soil texture fractions add up to {} but should '
                             'add up to approx 1'.format(total_fraction))

        self.bulk_den = bulk_den
        self.org_carb = org_carb
        self.f_clay = f_clay
        self.f_sand = f_sand
        self.f_silt = f_silt
        self.f_organic = f_organic
        self.soc_layer = soc_layer

        self.jules_soil_parameters = self.fill_parameters(
            f_sand, f_silt, f_clay, bulk_den, org_carb)

        self.jules_soil_parameters['hcon'] = \
            self._fill_thermal_conductivity_of_soil(
                self.f_sand, self.f_silt, self.f_clay, self['sm_sat'])

        self.jules_soil_parameters['hcap'] = \
            self._fill_heat_capacity_of_soil(
                self.f_sand, self.f_silt, self.f_clay, self['sm_sat'])

        if f_organic is not None:
            self.adapt_soil_parameters_for_organic_content(
                self.jules_soil_parameters, f_organic, soc_layer)

    def __getitem__(self, key):
        return self.jules_soil_parameters[key]

    @classmethod
    def from_soil_texture_class(cls, class_name, bulk_density, organic_carbon,
                                f_organic=None, soc_layer=None):
        """
        Generate soil type from the soil texture class.
        
        Parameters
        ----------
        soil_textural_class : str
            Soil textural class as in table 2 of [Cosby 1984]. See
            :obj:`SOIL_TEXTURAL_CLASS`.
         - bulk_den  : Bulk Density (kg / nm3)
         - org_carb  : Organic carbon (% weight)
         - f_organic : fraction of organic matter in soil
         - soc_layer : layer to take organic properties from (described in
           [Dankers 2011] table 2)

        """
        if class_name not in SOIL_TEXTURAL_CLASS:
            msg = 'Soil textural class sould be one of: {}'
            raise ValueError(msg.format(list(SOIL_TEXTURAL_CLASS.keys())))
        silt_percentage, sand_percentage, clay_percentage = (
            SOIL_TEXTURAL_CLASS[class_name])
        return cls(sand_percentage/100., silt_percentage/100.,
                   clay_percentage/100., bulk_density, organic_carbon,
                   f_organic=f_organic, soc_layer=soc_layer)

    @classmethod
    def _fill_thermal_conductivity_of_soil(cls, f_sand, f_silt, f_clay,
                                           sm_sat):
        '''
        Calculate thermal conductivity of non-organic soil (hcon)
        using [AncilDoc_SoilThermal] Method 1

        '''
        lambda_air  = 0.025  # W m^-1 K^-1
        lambda_clay = 1.16   # W m^-1 K^-1
        lambda_sand = 1.57   # W m^-1 K^-1
        lambda_silt = 1.57   # W m^-1 K^-1

        a = (1.0 - sm_sat) * f_clay
        b = (1.0 - sm_sat) * f_sand
        c = (1.0 - sm_sat) * f_silt

        return (
            lambda_air ** sm_sat
            * lambda_clay ** a * lambda_sand ** b * lambda_silt ** c)

    @classmethod
    def _fill_heat_capacity_of_soil(cls, f_sand, f_silt, f_clay, sm_sat):
        '''
        Calculate heat capacity of non-organic soil (hcap)
        using [AncilDoc_SoilThermal] Method 1
        
        '''
        c_clay = 2.373E6  # J m^-3 K^-1
        c_sand = 2.133E6  # J m^-3 K^-1
        c_silt = 2.133E6  # J m^-3 K^-1

        return (
            (1.0 - sm_sat)
            * (f_clay * c_clay + f_sand * c_sand + f_silt * c_silt))

    @classmethod
    def adapt_soil_parameters_for_organic_content(cls, soil_params, f_organic,
                                                  soc_layer):
        '''
        Adapts the soil parameters (calculated for non-organic)
        soil) to include organic content.
        '''
        if not (0.0 <= f_organic <= 1.0):
            raise ValueError('if given, f_organic should be 0.0 to 1.0')
        if soc_layer not in [1, 2, 3]:
            raise ValueError('soc_layer not recognised. Should be 1, 2 or 3.')

        i = soc_layer - 1

        # [Chadburn 2015] Equation A1
        soil_params['b'] = (
            (1.0 - f_organic) * soil_params['b'] +
            f_organic * ORGANIC_SOIL['b'][i]
            )
        # [Chadburn 2015] Equation A2
        soil_params['sathh'] = (
            soil_params['sathh'] ** (1.0 - f_organic) *
            ORGANIC_SOIL['sathh'][i] ** f_organic
            )
        # [Chadburn 2015] Equation A3
        soil_params['satcon'] = (
            soil_params['satcon'] ** (1.0 - f_organic) *
            ORGANIC_SOIL['satcon'][i] ** f_organic
            )
        # [Chadburn 2015] Equation A4
        soil_params['sm_sat'] = (
            (1.0 - f_organic) * soil_params['sm_sat'] +
            f_organic * ORGANIC_SOIL['sm_sat'][i]
            )
        # [Chadburn 2015] Equation A5
        # slight difference in psi to be consistent with the non-organic
        # calculation
        soil_params['sm_crit'] = cls._theta(
            0.033E6, soil_params['sm_sat'], soil_params['sathh'],
            soil_params['b'])

        # [Chadburn 2015] Equation A6
        # slight difference in psi to be consistent with the non-organic
        # calculation
        soil_params['sm_wilt'] = cls._theta(
            1.5E6, soil_params['sm_sat'], soil_params['sathh'],
            soil_params['b'])

        # [Chadburn 2015] Equation A7
        soil_params['hcap'] = (
            (1.0 - f_organic) * soil_params['hcap'] +
            f_organic * ORGANIC_SOIL['hcap'][i]
            )
        # [Chadburn 2015] Equation A8
        soil_params['hcon'] = (
            soil_params['hcon'] ** (1.0 - f_organic) *
            ORGANIC_SOIL['hcon'][i] ** f_organic
            )

    def __str__(self):
        '''
        prints soil parameters in a format that can be copy-and-pasted into the
        JULES_SOIL_PROPS namelist

        '''
        var = ''
        const_val = ''
        for key in self.PARAMETERS:
            var = var + ("'" + key + "'").rjust(14)
            const_val = const_val + '    ' + ("%10.4G" % (self[key]))
        return '{}\n{}'.format(var, const_val)


class VanGenuchten(_SoilType):
    @classmethod
    def _theta(cls, psi, sm_sat, sathh, b):
        '''Calculating vol. soil moisture.'''
        # h is the soil water pressure head
        h = psi / (RHO_W * G_GRAVITY)
        theta = 1.0 + (h / sathh)** ((1.0 + b) / b)
        return sm_sat * theta** (-1.0 / (1.0 + b))

    @classmethod
    def _psi(cls, theta, sm_sat, sathh, b):
        '''
        Calculates the matric water potential (psi) at a volumetric 
        soil moisture concentration (theta).
        abs_psi is in kg m^-1 s^-2 (i.e Pa). Positive.
        theta is m^3 water per m^3 soil
        '''
        h = -1.0 + (theta / sm_sat)** \
                         (-1.0 - b)
        h = sathh * h**(b / (1.0 + b))
        return h * (RHO_W * G_GRAVITY)

    @classmethod
    def fill_parameters(cls, f_sand, f_silt, f_clay, bulk_den, org_carb):
        '''
        Calculate soil parameters for non-organic soil
        using [AncilDoc_Cosby], which references [Cosby 1984]
        '''
        import rosetta_wrapper
        vg_res = rosetta_wrapper.get_van_genuchten_parameters(
            f_sand*100, f_silt*100, f_clay*100)
        res = {}

        # Residual soil moisture (cm3/cm3): Our modelling approach assumes that
        # we have subtracted this from the soil moisture state. So all you need
        # to do is subtract this value from the saturated soil moisture value
        # given by Rosetta.
        vg_res_soil_moist = vg_res[0]
        # Saturated soil moisture (cm3/cm3): We require (saturated - residual)
        # in units of m3/m3. However, as this is dimensionless, there is no
        # conversion factor required! So all you need to do is subtract the
        # residual soil moisture from this value and then set this as the
        # saturation value for JULES.
        vg_sat_soil_moist = vg_res[1]
        res['sm_sat'] = vg_sat_soil_moist - vg_res_soil_moist
        # Alpha (1/cm): We need the parameter 1/Alpha (m). So we need to
        # calculate 1.0 / (100.0 * Alpha). This value then replaces the
        # parameter 'sathh' in the ancil file. 
        vg_alpha = vg_res[2]
        res['sathh'] = 1.0 / (100.0 * vg_alpha)
        # n: This parameter is dimensionless, so there is no conversion factor.
        # However, we need the parameter 1.0 / (n - 1.0). This then replaces
        # the parameter 'b' in the ancil file
        vg_n = vg_res[3]
        res['b'] = 1.0 / (vg_n - 1.0)
        # Saturated hydraulic conductivity (cm/day). We need this parameter in
        # units of Kg/m2/s (equivalent to mm/s given the density of water is
        # 1000 Kg/m3). So you need to multiply the parameter by:
        # (10.0/86400.0).
        vg_sat_hydr_cond = vg_res[4]
        res['satcon'] = vg_sat_hydr_cond * (10.0/86400.0)

        # Volumetric soil moisture concentration at wilting point theta_w
        # This is calculated assuming a matric water potential (psi) of 1.5MPa
        # 1.0 MPa = 1.0E6 kg m^-1 s^-2
        # JULES_SOIL_PROPS units: m^3 water per m^3 soil
        res['sm_wilt'] = cls._theta(1.5E6, res['sm_sat'], res['sathh'],
                                    res['b'])

        # Volumetric soil moisture concentration at critical point theta_c
        # This is calculated assuming a matric water potential (psi)
        # of 0.033MPa.
        # JULES_SOIL_PROPS units: m^3 water per m^3 soil
        res['sm_crit'] = cls._theta(0.033E6, res['sm_sat'], res['sathh'],
                                    res['b'])
        
        # needed for l_triffid=T
        res['clay'] = f_clay
        
        # albsoil is not generated by this file - set to NaN
        res['albsoil'] = np.nan

        # Bulk density
        # Units in database are kg nm^-3, so need to multiply by
        # 1000.0 to convert to kg m^-3
        res['BD'] = 1000.0 * bulk_den

        # Soil carbon
        # Units in database are %, so need to convert to kg m^-3
        # (multiply by 1000.0 to get bulk density into kg m^-3 and
        # divide by 100.0 to get from percentage to fraction. 
        # Then assume all carbon is only in the topsoil layer, i.e., 0.3 m
        # (full calculation is equivalent to multiplying by 3.0)
        res['soil_carb'] = 3.0 * bulk_den * org_carb
        return res


class ClappHornberger(_SoilType):
    @classmethod
    def _theta(cls, psi, sm_sat, sathh, b):
        '''
        Calculates the volumetric soil moisture concentration (theta) at a
        matric water potential of psi using the brooks_and_corey_equation.
        psi is in kg m^-1 s^-2 (i.e Pa)
        theta is m^3 water per m^3 soil

        note: this differs from [AncilDoc_Cosby]
              because has extra factor sm_sat
              (checked this with Heather Ashton 1.4.15)
              
        See also [Best 2011] equation 58.      
        '''

        # h is the soil water pressure head
        h = psi / (RHO_W * G_GRAVITY)
        return sm_sat * ((sathh / h) ** (1.0 / b))

    @classmethod
    def _psi(cls, theta, sm_sat, sathh, b):
        '''
        Calculates the matric water potential (psi) at a volumetric 
        soil moisture concentration (theta).  Inverse
        brooks_and_corey_equation.
        abs_psi is in kg m^-1 s^-2 (i.e Pa). Positive.
        theta is m^3 water per m^3 soil    

        '''
        return sathh * (theta / sm_sat) ** (- b) * (RHO_W * G_GRAVITY)

    @classmethod
    def fill_parameters(cls, f_sand, f_silt, f_clay, bulk_den, org_carb):
        '''
        Calculate soil parameters for non-organic soil
        using [AncilDoc_Cosby], which references [Cosby 1984]
        '''
        res = {}

        # Clapp Hornberger parameter b
        # JULES_SOIL_PROPS units: dimensionless.
        res['b'] = 3.10 + 15.70 * f_clay - 0.3 * f_sand

        # Saturated soil water suction SATHH in terms of log to the base 10
        # (n.b. used to use nat log)
        # JULES_SOIL_PROPS units: m. (says 'the *absolute* value of the soil
        # matric suction at saturation')
        res['sathh'] = 0.01 * (10.0 ** (2.17 - 0.63 * f_clay - 1.58 * f_sand))

        # Saturated hydrological soil conductivity Ks in terms of log to the
        # base 10 (n.b. used to use nat log)
        # JULES_SOIL_PROPS units: kg m^-2 s^-1
        res['satcon'] = 10.0 ** (-2.75 - 0.64 * f_clay + 1.26 * f_sand)

        # Volumetric soil water concentration at saturation point theta_s
        # JULES_SOIL_PROPS units: m^3 water per m^3 soil
        res['sm_sat'] = 0.505 - 0.037 * f_clay - 0.142 * f_sand

        # Volumetric soil moisture concentration at wilting point theta_w
        # This is calculated assuming a matric water potential (psi) of 1.5MPa
        # 1.0 MPa = 1.0E6 kg m^-1 s^-2
        # JULES_SOIL_PROPS units: m^3 water per m^3 soil
        res['sm_wilt'] = cls._theta(1.5E6, res['sm_sat'], res['sathh'],
                                    res['b'])

        # Volumetric soil moisture concentration at critical point theta_c
        # This is calculated assuming a matric water potential (psi)
        # of 0.033MPa.
        # JULES_SOIL_PROPS units: m^3 water per m^3 soil
        res['sm_crit'] = cls._theta(0.033E6, res['sm_sat'], res['sathh'],
                                    res['b'])
        
        # needed for l_triffid=T
        res['clay'] = f_clay
        
        # albsoil is not generated by this file - set to NaN
        res['albsoil'] = np.nan

        # Bulk density
        # Units in database are kg nm^-3, so need to multiply by
        # 1000.0 to convert to kg m^-3
        res['BD'] = 1000.0 * bulk_den

        # Soil carbon
        # Units in database are %, so need to convert to kg m^-3
        # (multiply by 1000.0 to get bulk density into kg m^-3 and
        # divide by 100.0 to get from percentage to fraction. 
        # Then assume all carbon is only in the topsoil layer, i.e., 0.3 m
        # (full calculation is equivalent to multiplying by 3.0)
        res['soil_carb'] = 3.0 * bulk_den * org_carb
        return res


if __name__ == '__main__':
    f_clay, f_sand, f_silt, bulk_den, org_carb = 0.2, 0.3, 0.5, 1.5, 2 
    print(VanGenuchten(f_clay, f_sand, f_silt, bulk_den, org_carb))

    print(ClappHornberger(f_clay, f_sand, f_silt, bulk_den, org_carb))
