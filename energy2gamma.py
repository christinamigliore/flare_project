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

# Electron mass in kilograms
E_MASS_KG = 9.1094*10**-31
# Proton mass in kilograms
P_MASS_KG = 1.6726*10**-27
# Speed of light [m/s]
SPEED_OF_LIGHT = 2.9979*10**8
# converts ergs to Joules
ERGS_2_JOULES = 10**-7
# converts keV to Joules
KEV_2_JOULES = 1.6022*10**-16

def converting_energies(energy_type, energy_value):
    """
    Converts energies from kev, mev, or ergs to joules
    """
    if energy_type == 'ergs':
        energy_value = energy_value*ERGS_2_JOULES
    if energy_type == 'joules':
        energy_value = energy_value
    if energy_type == 'mev':
        energy_value == energy_value*10**3*KEV_2_JOULES
    else:
        energy_value = energy_value*KEV_2_JOULES
    return energy_value

def finding_gamma(energy_value, energy_type='kev', particle_type='electron', total_energy=True):
    """
    Calculates the lorentz gamma using inputed energy
    """
    energy_value = converting_energies(energy_type, energy_value)
    if particle_type == 'proton':
        mass = P_MASS_KG
    else:
        mass = E_MASS_KG

    if total_energy:
        lorentz_gamma = energy_value/(mass*SPEED_OF_LIGHT**2)
    else:
        lorentz_gamma = energy_value/(mass*SPEED_OF_LIGHT**2) + 1
    return lorentz_gamma
