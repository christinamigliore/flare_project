"""
PURPOSE
-------
	Calculate the b limit value in particle-particle
    interaction

INPUTS
------
	energy_value: float
		Particle (assumed to be an electron) energy
		(assumed to be in keV and assumed to be the kinetic energy)

	particle_type: string, optional 
		Used to finding the mass of the particle, can be proton or electron
		(If not specified, the mass is set to the mass of an electron)

OUTPUTS
-------
	b_limit: float

"""
import unittest

import energy2vel as e

# Electron mass in grams
E_MASS_G = 9.1094*10**-28
# Proton mass in grams
P_MASS_G = 1.6726*10**-24
# charge of electron in coulombs
E_CHARGE = 1.602*10**-19


psi_min = 1*10**-4

def get_b_limit(energy_value, particle_type='proton'):
    """
    Calculate the b limit value in particle-paritcle interactions
    """
    if energy_value == 0:
        energy_value == 15

    v = e.energy2vel(energy_value)
    m_test = E_MASS_G
    if particle_type == 'proton':
        m_target = P_MASS_G
    else:
        m_target = E_MASS_G

    mu_rm = (m_test*m_target)/(m_test+m_target)

    b_limit = (2*E_CHARGE*E_CHARGE)/(mu_rm*v*v*psi_min)
    return b_limit

class Test_get_blimit(unittest.TestCase):
    def test_get_b_limit(self):
        self.assertEqual(get_b_limit(600),  1.8004325017426443e-27)

if __name__ == '__main__':
    unittest.main()
