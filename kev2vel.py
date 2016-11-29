"""
PURPOSE
-------
	Calculate an electron's or particle's speed based on
	its kinetic energy.

INPUTS
------
	energy_value: float
		Particle (assumed to be an electron) energy
		(assumed to be in keV and assumed to be the kinetic energy)

	particle_type: string
		Used to finding the mass of the particle, can be proton or electron

OUTPUTS
-------
	v_total: float
		speed in (cm/s)

"""
import math

# CONSTANTS

#1 keV=1.602e-9 ergs
# converts keV to ergs
KEV_2_ERGS = 1.6022*10**-9
# Electron mass in grams
E_MASS_G = 9.1094*10**-28
# Proton mass in grams
P_MASS_G = 1.6726*10**-24

def kev_2_vel(energy_value, particle_type='electron'):
    """
    Calculates total velocity from inputed kinetic energy (keV)
    """
    if particle_type == 'proton':
        mass = P_MASS_G
    else:
        mass = E_MASS_G
    v_total = math.sqrt((2*energy_value*KEV_2_ERGS)/mass)
    return v_total