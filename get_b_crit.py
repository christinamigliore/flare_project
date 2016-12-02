"""
PURPOSE
-------
	Calculate the b critical value in particle-particle
    interaction

INPUTS
------
	energy_value: float
		Particle (assumed to be an electron) energy
		(assumed to be in keV and assumed to be the kinetic energy)

OUTPUTS
-------
	b_limit: float

"""
#1 keV=1.602e-9 ergs
# converts keV to ergs
KEV_2_ERGS = 1.6022*10**-9
# Electron mass in grams
E_MASS_G = 9.1094*10**-28
# Proton mass in grams
P_MASS_G = 1.6726*10**-24
# charge of electron in coulombs
E_CHARGE = 1.602*10**-19

def get_b_crit(energy_value):
    """
    Calculates the b critical value
    """
    if energy_value ==  0:
        energy_value = 15
    
    b_crit = (E_CHARGE*E_CHARGE)/(energy_value*KEV_2_ERGS)
    return b_crit