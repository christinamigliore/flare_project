import math
import numpy as np
from astropy import units as u
from scipy import integrate

# Charge of electron in esmu
E_CHARGE = 4.80e-10
# Electron mass in grams
E_MASS_G = 9.1094e-28
# Proton mass in grams
P_MASS_G = 1.6726219e-24
# Electron density in cm^3
E_DENSITY = 1e9
# erg^2 to keV^2
conversion_factor = 3.896e17

# FIX THIS, EQ ONLY VALID FOR T > 1E5 K
def calculate_coulomb_log(electron_n=1e9*u.cm**-3, Te=2e6*u.K):
    electron_n = electron_n.to(u.cm**-3)
    if Te.unit == 'eV':
        Te = Te.to(units.K, equivalencies=units.temperature_energy())
    coulomb_log = np.log((8e6*electron_ke.value)/np.sqrt(electron_n.value))
    return coulomb_log

def collisions(energy_value, v_int, electron_n=1e9*u.cm**-3, Te=2e6*u.K):
    """
    PURPOSE
    -------
        Function calculates the distance "N_stop" that a particle's energy turns 
        completely into thermal energy.

    INPUTS
    ------
        energy_value: float
            Particle (assumed to be an electron) energy
            (assumed to be in keV and assumed to be the kinetic energy)

        v_init: array
            The initial velocity of particle in component form. [Vx0, Vy0, Vz0]

    OUTPUTS
    -------
        N_stop: float
            Distance where all of the particle's energy has turned to thermal energy
    """
    if energy_value.value ==  0:
        energy_value = 15 * u.keV

    energy_value = energy_value.to(u.keV)

    vx = v_int[0] * u.m/u.s
    vy = v_int[1] * u.m/u.s
    vz = v_int[2] * u.m/u.s

    vx = vx.to(u.cm/u.s)
    vy = vy.to(u.cm/u.s)
    vz = vz.to(u.cm/u.s)

    velocity_mag = np.sqrt(vx.value**2 + vy.value**2 + vz.value**2)
    mu_naught = vz.value/velocity_mag
    coulomb_log = calculate_coulomb_log(electron_n=electron_n, Te=Te)
    C = 2*np.pi*(E_CHARGE**4)*coulomb_log
    N_stop = (mu_naught*energy_value.value**2)/(3*C*conversion_factor)
    return N_stop

def change_in_energy(X, n, mu_naught, C):
    """
    PURPOSE
    -------
        Function calculates the change in energy particle due to collisions 
        in thermal plasma.

    INPUTS
    ------
        energy_value: float
            Particle (assumed to be an electron) energy
            (assumed to be in keV and assumed to be the kinetic energy)

        constant: float
            The value found in the function "collisons" which corresponds to 
            the plasma parameters

        stop_distance: float
            The distance at which to stop calculating the change in energy.

    OUTPUTS
    -------
        de_dn: float
            Change in energy for that distance traveled.
    """
    de_dn = (-1*C)/(mu_naught*X*(1 - (3*C*n)/(mu_naught*X**2))**float(2/3)) 

    return de_dn


def calculate_change_in_energy(energy_value, v_int, dn, electron_n=1e9*u.cm**-3, Te=2e6*u.K):
    """
    PURPOSE
    -------
        Function calculates the change in energy at each distance step

    INPUTS
    ------
        energy_value: float
            Particle's (assumed to be an electron) initial energy
            (assumed to be in keV and assumed to be the kinetic energy)


        v_init: array
            The initial velocity of particle in component form. [Vx0, Vy0, Vz0]

        dn: float
            The distance step  

    OUTPUTS
    -------
        n: list
            List of the energies of the particle at each distance step
    """
    vx = v_int[0] * u.m/u.s
    vy = v_int[1] * u.m/u.s
    vz = v_int[2] * u.m/u.s

    vx = vx.to(u.cm/u.s)
    vy = vy.to(u.cm/u.s)
    vz = vz.to(u.cm/u.s)

    velocity_mag = np.sqrt(vx.value**2 + vy.value**2 + vz.value**2)
    energy_value = energy_value.to(u.erg)
    mu_naught = vz.value/velocity_mag
    coulomb_log = calculate_coulomb_log(electron_n=electron_n, Te=Te)
    C = 2*np.pi*(E_CHARGE**4)*coulomb_log

    stop_distance = collisions(energy_value, v_int)
    ns = np.arange(0, 1e10, dn)
    new_energy = integrate.odeint(change_in_energy, energy_value.value, ns, args=(mu_naught, C))
    new_energy = new_energy * u.erg
    new_energy = new_energy.to(u.keV)

    return new_energy
