"""
Common thermodynamics.

1. Saturation vapour pressure.

2. Dry adiabat.

3. Parcel ascent.

"""

from constants import SPECIFIC_GAS_CONSTANT_FOR_DRY_AIR as RD
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
            225–233, doi:10.1175/BAMS-86-2-225.

    """
    # Convert to degrees C
    ta = ta - 273.15
    c1 = 610.94 # Pa
    a1 = 17.625
    b1 = 243.04 # C
    return c1 * np.exp(a1*t/(b1+t))
