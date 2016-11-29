"""
PURPOSE
-------
	Calculate an electron's or particle's speed based on
	its kinetic energy.  Relativistic effects are taken
	into account for gamma =>1.2

INPUTS
------
	energy_value: float
		Particle (assumed to be an electron) energy
		(assumed to be in keV and assumed to be the kinetic energy)

	total_energy: boolean, optional
		True if the input was the total energy of the particle
        not just the kinetic energy
        (If not specified, set to True)

	energy_type: string, optional
		Options: kev, joules, ergs
		Converts energy into different units
		(If not specified, the energy type is set to kev)

	particle_type: string, optional
		Used to finding the mass of the particle, can be proton or electron
		(If not specified, the mass is set to the mass of an electron)

OUTPUTS
-------
	v_total: float
		speed in (m/s)

"""
import math

import energy2gamma as e

# CONSTANTS

# Electron mass in kilograms
E_MASS_KG = 9.1094*10**-31
# Proton mass in kilograms
P_MASS_KG = 1.6726*10**-27
# Speed of light [m/s]
SPEED_OF_LIGHT = 2.9979*10**8


def energy2vel(energy_value, energy_type='kev', particle_type='electron', total_energy=True):
    """
    Calculates the velocity from the inputed energy
    """
    if particle_type == 'proton':
        mass = P_MASS_KG
    else:
        mass = E_MASS_KG

    rest_energy = mass*SPEED_OF_LIGHT**2
    lorentz_gamma = e.finding_gamma(energy_value, energy_type, particle_type, total_energy)
    energy_value = e.converting_energies(energy_type, energy_value)

    if lorentz_gamma < 1.2:
        if total_energy:
            v_total = math.sqrt(2*(energy_value - rest_energy)/mass)
        else:
            v_total = math.sqrt((2*energy_value)/mass)
    else:
        if total_energy:
            v_total = SPEED_OF_LIGHT*((1-(((energy_value-rest_energy)
                                           /(rest_energy)+1)**(-2)))**(0.5))
        else:
            v_total = SPEED_OF_LIGHT*((1-((energy_value/(rest_energy)+1)**(-2)))**(0.5))
    return v_total
