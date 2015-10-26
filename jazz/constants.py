"""Thermodynamic constants, largely taken from Wallace and Hobbs Ch 3."""

AVOGRADROS_NUMBER = 6.022E23
UNIVERSAL_GAS_CONSTANT = 8.3145  # J K-1 mol-1 


MOLECULAR_WEIGHT_FOR_DRY_AIR = 28.97  # g
MOLECULAR_WEIGHT_FOR_WATER = 18.016   # g
SPECIFIC_GAS_CONSTANT_FOR_DRY_AIR = (UNIVERSAL_GAS_CONSTANT /
                                     MOLECULAR_WEIGHT_FOR_DRY_AIR)
SPECIFIC_GAS_CONSTANT_FOR_WATER_VAPOUR = (UNIVERSAL_GAS_CONSTANT /
                                          MOLECULAR_WEIGHT_FOR_WATER)

GRAVITATIONAL_ACCELERATION = 9.81 # m s-2

# Molecular weights
MASS_AIR = 28.97
MASS_CO2 = 44.01
MASS_CH4 = 16.04
MASS_N2O = 44.01
MASS_HFC125 = 102.02
MASS_HFC134A = 102.03
MASS_CFC11 = 137.37
MASS_CFC12 = 120.91
MASS_CFC113 = 187.38
MASS_CFC114 = 170.92
MASS_HCFC22 = 86.47
