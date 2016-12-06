"""
PURPOSE
-------
	Calculatek2 an electron's or particle's speed based on
	its kinetic energy.

INPUTS
------
	energy_value: float
		Particle (assumed to be an electron) energy
		(I6n keV and assumed to be the kinetic energy)

	particle_type: string
		Used to finding the mass of the particle, can be proton or electron

OUTPUTS
-------
	v_total: float
		speed in (cm/s)

"""
import math
import unittest

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

class Test_kev_2_vel(unittest.TestCase):
    def test_kev2vel(self):
        self.assertEqual(kev_2_vel(250), 29655037635.500153)
        self.assertEqual(kev_2_vel(30, 'proton') 239738587.75334778)

if __name__ == '__main__':
    unittest.main()
