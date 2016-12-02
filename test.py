"""
PURPOSE
-------
	Unit test each of the five programs built to 
	calculate b_limit, b_critical, velcity, lorentz gamma

INPUTS
------
	None

OUTPUTS
-------
	"OK" if tests have passed
	"ERROR" if there is an error in the code
	"FAILED" if one or more of the tests failed
"""
import unittest

import energy2gamma as e 
import energy2vel as v
import kev2vel as vel
import get_b_crit as g
import get_b_limit as b
 
class Test_energy_2_gamma(unittest.TestCase):
 
    def test_converting_energies(self):
        self.assertEqual(e.converting_energies('kev', 600), 9.613200000000001e-7)
        self.assertEqual(e.converting_energies('ergs', 9.6132e-7), 9.6132e-7)
        self.assertEqual(e.converting_energies('joules', 9.61306e-14), 9.613059999999999e-21)
        self.assertEqual(e.converting_energies('mev', 0.6), 9.613200000000001e-7)

    def test_finding_gamma(self):
        self.assertEqual(e.finding_gamma(600), 1.1742049878397725)
        self.assertEqual(e.finding_gamma(600, 'kev', 'proton'), 0.0006395015494575886)
        self.assertEqual(e.finding_gamma(600, 'kev', 'electron', False) 2.1742049878397722)
        self.assertEqual(e.finding_gamma(600, 'kev', 'proton', False), 1.0006395015494576)

    def test_energy2vel(self):
        self.assertEqual(v.energy2vel(600), 17695483468.69996)
        self.assertEqual(v.energy2vel(700), 20489792344.582314)
        self.assertEqual(v.energy2vel(700, 'kev', 'electron', False), 27179426697.076885)
        self.assertEqual(v.energy2vel(600, 'kev', 'electron', False), 26619880387.886257)
        self.assertEqual(v.energy2vel(600, 'kev', 'proton', False), 1072143558.0925685)
        self.assertEqual(v.energy2vel(700, 'kev', 'proton', False), 1158047398.5777202)

    def test_kev2vel(self):
        self.assertEqual(vel.kev_2_vel(250), 29655037635.500153)
        self.assertEqual(vel.kev_2_vel(30, 'proton') 239738587.75334778)

    def test_get_b_crit(self):
        self.assertEqual(g.get_b_crit(600), 2.66966670827612e-32)

    def test_get_b_limit(self):
        self.assertEqual(b.get_b_limit(600),  1.8004325017426443e-27)
 
if __name__ == '__main__':
    unittest.main()
