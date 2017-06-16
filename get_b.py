import unittest

import energy2vel as e

#1 keV=1.602e-9 ergs
# converts keV to ergs
KEV_2_ERGS = 1.6022*10**-9
# Electron mass in grams
E_MASS_G = 9.1094*10**-28
# Proton mass in grams
P_MASS_G = 1.6726*10**-24
# charge of electron in coulombs
E_CHARGE = 1.602*10**-19
psi_min = 10**-4

def get_b_crit(energy_value):
    """
    PURPOSE
    -------
        Calculates the b critical value in particle-particle
        interaction.

    INPUTS
    ------
        energy_value: float
            Particle's (assumed to be an electron) energy
            (assumed to be in keV and assumed to be the kinetic energy).

    OUTPUTS
    -------
        b_limit: float
    """
    if energy_value ==  0:
        energy_value = 15
    
    b_crit = (E_CHARGE*E_CHARGE)/(energy_value*KEV_2_ERGS)
    return b_crit

def get_b_limit(energy_value, particle_type='proton'):
    """
    PURPOSE
    -------
        Calculate the b limit value in particle-particle interaction.

    INPUTS
    ------
        energy_value: float
            Particle (assumed to be an electron) energy
            (assumed to be in keV and assumed to be the kinetic energy)

        particle_type: string (optional) 
            Used to finding the mass of the particle, can be proton or electron
            (If not specified, the mass is set to the mass of an electron)

    OUTPUTS
    -------
        b_limit: float
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

    b_limit = (2*E_CHARGE*E_CHARGE)/(mu_rm*psi_min*v**2)
    return b_limit

class Test_energy_2_gamma(unittest.TestCase):
    def test_get_b_crit(self):
        self.assertEqual(get_b_crit(600), 2.66966670827612e-32)

class Test_get_blimit(unittest.TestCase):
    def test_get_b_limit(self):
        self.assertEqual(get_b_limit(600),  1.8004325017426443e-27)

if __name__ == '__main__':
    unittest.main()