"""
The Debye length is the distance to which a point charge can have its influence felt.
"""

import math

def calculate_debye(electron_nd=10**32, ion_nd=10**32, electron_ke=2*10**6,
                    ion_ke=2*10**6):
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
            The Debye length (meters)

    """
    k_boltzmann = 1.38064852 * 10**-23
    electron_charge = 1.6021766208 * 10**-19
    total_ke = (electron_nd*electron_ke + ion_nd*ion_ke)/(electron_nd + ion_nd)
    debye = math.sqrt(k_boltzmann*total_ke/(4*math.pi*electron_nd*electron_charge**2))
    return debye

if __name__ == '__main__':
    print(calculate_debye())
