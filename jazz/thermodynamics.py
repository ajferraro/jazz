"""
Common thermodynamics.

1. Saturation vapour pressure.

2. Dry adiabat.

3. Parcel ascent.

"""

import numpy as np

import cf_units
import iris

from constants import SPECIFIC_HEAT_CONSTANT_PRESSURE as CP
from constants import LATENT_HEAT_OF_CONDENSATION_OF_WATER as LV
from constants import SPECIFIC_GAS_CONSTANT_FOR_DRY_AIR as RD
from constants import SPECIFIC_GAS_CONSTANT_FOR_WATER_VAPOUR as RV
from constants import EPSILON
from constants import GRAVITATIONAL_ACCELERATION as GRAV


def parcel_ascent():
    pass


def saturated_adiabat():
    pass


def dry_adiabat():
    pass


def svp(ta):
    """Calculate saturation vapour pressure using the
    August-Roche-Magnus formula (Lawrence 2005).

    Args:
        ta: temperature in Kelvin

    References:
        Lawrence, M. G., 2005: The relationship between relative
            humidity and the dewpoint temperature in moist air: A simple
            conversion and applications. Bull. Am. Meteorol. Soc., 86,
            225-233, doi:10.1175/BAMS-86-2-225.

    """
    if isinstance(ta, iris.cube.Cube):
        is_cube = True
        ta_data = ta.data
    else:
        is_cube = False
        ta_data = ta
    # Convert to degrees C
    ta_data = ta_data - 273.15
    c1 = 610.94 # Pa
    a1 = 17.625
    b1 = 243.04 # C
    svp_data = c1 * np.exp(a1*ta_data/(b1+ta_data))
    if is_cube:
        svp_cube = ta.copy(data=svp_data)
        svp_cube.units = cf_units.Unit('Pa')
        svp_cube.rename('water_vapor_partial_pressure_in_air')
        return svp_cube
    else:
        return svp_data


def hus_from_hur(hur, ta, pres):
    """Calculate specific humidity from relative humidity. Uses equation 3.64
    of Wallace and Hobbs (2006).

    Args:
        hur (float): relative humidity.
        ta (float): air temperature.
        pres (float): air pressure.

    Returns:
        float

    """
    # Calculate saturation humidity mixing ratio
    sat_vapour_pres = svp(ta)
    sat_humidity_mixing_ratio = shmr(sat_vapour_pres, pres)

    # Calculate humidity mixing ratio
    humidity_mixing_ratio = sat_humidity_mixing_ratio*hur/100
    print 'Converting hur from % to fraction'

    # Convert to specific humidity
    return hus_from_mr(humidity_mixing_ratio)


def hur_from_hus(hus, ta, pres):
    """Calculate relative humidity from specific humidity.  Uses Equation 3.64
    of Wallace and Hobbs (2006).

    Args:
        hus (float): specific humidity.
        ta (float): air temperature.
        pres (float): air pressure.

    Returns:
        float

    """
    # Convert specific humidity into humidity mixing ratio
    mr = mr_from_hus(hus)

    # Calculate sat humidity mixing ratio
    sat_vapour_pres = svp(ta)
    mr_sat = shmr(sat_vapour_pres, pres)

    print 'Converting hur from fraction to %'
    return 100*mr/mr_sat


def hus_from_mr(mr):
    """Calculate specific humidity from mixing ratio."""
    return mr/(1+mr)


def mr_from_hus(hus):
    """Calculate mixing ratio from specific humidity."""
    return hus/(1+(-1*hus))


def vp(mr, pres):
    """Calculate water vapour pressure from Equation 3.59 of Wallace and Hobbs
    (2008).

    Args:
        mr (float): mixing ratio.
        pres (float): air pressure.

    Returns:
        float

    """
    return (mr*pres)/(mr+EPSILON)


def virtual_temperature(ta, mr):
    """Calculate virtual temperature using a more exact form of Equation 3.60
    of Wallace and Hobbs (2008).

    Args:
        ta (float): air temperature.
        mr (float): mixing ratio.

    Returns:
        float

    """
    return ta*(mr+EPSILON)/(EPSILON*(1+mr))


def calc_shmr(sat_vapour_pres, pres):
    """Calculate saturation humidity mixing ratio using Equation 3.63 of
    Wallace and Hobbs (2008).

    Args:
        sat_vapour_pres (float): saturation vapour pressure.
        pres (float): air pressure.

    Returns:
        float

    """
    return EPSILON*(sat_vapour_pres/
                    (pres-sat_vapour_pres))


def potential_temperature(ta, pres, pres_ref=1E5):
    """Calculate potential temperature.

    Args:
        ta (float): air temperature (K).
        pres (float): air pressure (hPa) corresponding to `ta`.
        pref_ref (Optional[float]): reference pressure (hPa).

    Returns:
        float

    """
    return ta*(pres_ref/pres)**(RD/CP)


def calc_dewpoint(ta, hurs):
    A = 17.625
    B = iris.cube.Cube(243.04, units='degrees_Celsius')

    ta_celsius = ta.copy()
    ta_celsius.convert_units('degrees_Celsius')
    print ta_celsius.units
    loghurs = hurs.copy(data=np.log(hurs.data))
    numerator = (loghurs + (A*ta_celsius)/(B.data+ta_celsius)) * B
    denominator = A + -1*loghurs + -1*(A*ta_celsius)/(B.data+ta_celsius)
    print numerator.units
    print denominator.units
    td = numerator/denominator
    td.convert_units('Kelvin')
    return td


def lcl(ts, ps, hurs):
    mrsat = calc_shmr(svp(ts), ps)
    mr = hurs*mrsat
    td = calc_dewpoint()
