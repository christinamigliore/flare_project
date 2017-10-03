"""
The Debye length is the distance to which a point charge can have its influence felt.
"""

import math

def calculate_debye(electron_nd=1e15, ion_nd=1e15, electron_ke=2e6,
                    ion_ke=2e6):
    """
    Calculates the Debye length using the formula found in "The Physics of
    Solar Flares" from chapter 3, page 61.

    Parameters
    ----------
	    electron_nd : float (optional)
	        Electron number density

		electron_ke : float (optional)
		    Electron kinetic energy

		ion_nd : float (optional)
		    Ion number density

		ion_ke : float (optional)
		    Ion kinetic energy

    Returns
    -------
        debye : float
            The Debye length (m)

    """
    K_BOLTZMANN = 1.3806e-23
    ELECTRON_CHARGE = 1.6021766208e-19
    EPSILON = 8.8541e-12
    total_ke = (electron_nd*electron_ke + ion_nd*ion_ke)/(electron_nd + ion_nd)
    print "test"
    debye = math.sqrt(EPSILON*K_BOLTZMANN*total_ke/(electron_nd*(ELECTRON_CHARGE**2)))
    return debye

if __name__ == '__main__':
    print(calculate_debye())