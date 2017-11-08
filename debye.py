"""
The Debye length is the distance to which a point charge can have its influence felt.
Written in SI units.
"""

import math
import numpy as np
from astropy import units as u

K_BOLTZMANN = 1.3806e-23
ELECTRON_CHARGE = 1.6021766208e-19
EPSILON = 8.8541e-12

def calculate_debye(electron_nd=1e15*u.m**-3, ion_nd=1e15*u.m**-3, Te=2e6*u.K,
                    Ti=2e6*u.K):
    """
    Calculates the Debye length using the formula found in "The Physics of
    Solar Flares" from chapter 3, page 61.

    Parameters
    ----------
	    electron_nd : float (optional)
	        Electron number density (m**-3)

		Te : float (optional)
		    Electron temperature (K)

		ion_nd : float (optional)
		    Ion number density (m**-3)

		Ti : float (optional)
		    Ion temperature (K)

    Returns
    -------
        debye : float
            The Debye length (m)
    """
    electron_nd = electron_nd.to(u.m**-3)
    ion_nd = ion_nd.to(u.m**-3)

    if Te.unit == 'eV':
        Te = Te.to(units.K, equivalencies=units.temperature_energy())
    if Ti.unit == 'ev':
        Ti = Ti.to(units.K, equivalencies=units.temperature_energy())

    total_ke = (electron_nd.value*Te.value + ion_nd.value*Ti.value)/(electron_nd.value + ion_nd.value)
    debye = math.sqrt(EPSILON*K_BOLTZMANN*total_ke/(electron_nd.value*ELECTRON_CHARGE**2))
    return debye * u.m

if __name__ == '__main__':
    print(calculate_debye())
