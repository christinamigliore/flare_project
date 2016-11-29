"""
PURPOSE
-------
	Given the energy of a particle, determine the Lorentz gamma

INPUTS
------
	energy_value: float
		Particle (assumed to be an electron) energy
		(assumed to be in keV and assumed to be the kinetic energy)

	energy_type: string, optional
		Options: kev, mev joules, ergs
		Converts energy into joules
		(If not specified, the energy type is set to kev)

	particle_type: string, optional
		Used to finding the mass of the particle, can be proton or electron
		(If not specified, the mass is set to the mass of an electron)

	total_energy: boolean, optional
		True if the input was the total energy of the particle
        not just the kinetic energy
        (If not specified, set to True)


OUTPUTS
-------
	lorentz_gamma: float
		Dimensionless Lorentz gamma
"""

# CONSTANTS

# Electron mass in grams
E_MASS_G = 9.1094*10**-28
# Proton mass in grams
P_MASS_G = 1.6726*10**-24
# Speed of light [cm/s]
SPEED_OF_LIGHT = 2.9979*10**10
# converts keV to ergs
KEV_2_ERGS = 1.6022*10**-9

def converting_energies(energy_type, energy_value):
    """
    Converts energy value from kev, mev, or joules to ergs
    """
    if energy_type == 'ergs':
        energy_value = energy_value
    if energy_type == 'joules':
        energy_value = energy_value*10**-7
    if energy_type == 'mev':
        energy_value == energy_value*10**3*KEV_2_ERGS
    else:
        energy_value = energy_value*KEV_2_ERGS
    return energy_value

def finding_gamma(energy_value, energy_type='kev', particle_type='electron', total_energy=True):
    """
    Calculates the lorentz gamma using inputed energy
    """
    energy_value = converting_energies(energy_type, energy_value)
    if particle_type == 'proton':
        mass = P_MASS_G
    else:
        mass = E_MASS_G

    if total_energy:
        lorentz_gamma = energy_value/(mass*SPEED_OF_LIGHT**2)
    else:
        lorentz_gamma = energy_value/(mass*SPEED_OF_LIGHT**2) + 1
    return lorentz_gamma
