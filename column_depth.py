import math
import numpy as np
from astropy import units as u
from scipy import integrate
import debye

# Charge of electron in esmu
E_CHARGE_ESMU = 4.80e-10
# Charge of electron in SI
E_CHARGE_SI = 1.6021e-19*u.C
# Electron mass in grams
E_MASS_G = 9.1094e-28*u.g
# Proton mass in grams
P_MASS_G = 1.6726219e-24*u.g
#Boltzmann Constant in SI
k_B = 1.38064852e-23*((u.m**2*u.kg)/(u.K*u.s**2))
#Vacuum of permeability
eps0 = 8.854e-12

def calculate_coulomb_log(electron_nd=1e9*u.cm**-3, ion_nd=1e9*u.cm**-3, Te=2e6*u.K, Ti=2e6*u.K):
    """
    PURPOSE
    -------
        Function calculates the Coulomb Logrithm through the plasma parameters.
        Note: input's units will automatically be converted into SI (need to be in 
        the astropy.units format).

    INPUTS
    ------
        electron_nd : float (optional)
            Electron number density

        ion_nd : float (optional)
            Ion number density

        Te : float (optional)
            Electron temperature

        Ti : float (optional)
            Ion temperature

        Returns
    -------
        coulomb_log : float
            The Coulomb logrithm (dimensionless)
    """
    electron_nd = electron_nd.to(u.m**-3)
    ion_nd = ion_nd.to(u.m**-3)
    if Te.unit == 'eV':
        Te = Te.to(units.K, equivalencies=units.temperature_energy())
    if Ti.unit == 'eV':
        Ti = Ti.to(units.K, equivalencies=units.temperature_energy())
    E_MASS = E_MASS_G.to(u.kg)
    P_MASS = P_MASS_G.to(u.kg)

    b_max = debye.calculate_debye(electron_nd=electron_nd, ion_nd=ion_nd, Te=Te, Ti=Ti)
    reduced_mass = E_MASS * P_MASS / (E_MASS + P_MASS)
    V = np.sqrt(3 * k_B * Te / reduced_mass)
    b_min = E_CHARGE_SI * E_CHARGE_SI / (4 * np.pi * eps0 * reduced_mass * V**2)
    coulomb_log = np.log(b_max.value/b_min.value)

    return coulomb_log

def collisions(energy_value, v_int, electron_nd=1e9*u.cm**-3, ion_nd=1e9*u.cm**-3, Te=2e6*u.K, Ti=2e6*u.K):
    """
    PURPOSE
    -------
        Function calculates the distance "N_stop" where the particle's energy turns 
        completely into thermal energy.

    INPUTS
    ------
        energy_value: float
            Particle (assumed to be an electron) energy
            (assumed to be in keV and assumed to be the kinetic energy)

        v_init: array
            The initial velocity of particle in component form. [Vx0, Vy0, Vz0]

        electron_nd : float (optional)
            Electron number density

        ion_nd : float (optional)
            Ion number density

        Te : float (optional)
            Electron temperature

        Ti : float (optional)
            Ion temperature

    OUTPUTS
    -------
        N_stop: float
            Distance where all of the particle's energy has turned to thermal energy
    """
    if energy_value.value ==  0:
        energy_value = 15 * u.keV

    energy_value = energy_value.to(u.erg)

    vx = v_int[0].to(u.cm/u.s)
    vy = v_int[1].to(u.cm/u.s)
    vz = v_int[2].to(u.cm/u.s)

    velocity_mag = np.sqrt(vx**2 + vy**2 + vz**2)
    mu_naught = vz/velocity_mag
    coulomb_log = calculate_coulomb_log(electron_nd=electron_nd, ion_nd=ion_nd, Te=Te, Ti=Ti)
    C = 2*np.pi*(E_CHARGE_ESMU**4)*coulomb_log
    N_stop = (mu_naught*(energy_value.value**2))/(3*C)

    return N_stop.value*u.cm

def change_in_energy(X, n, mu_naught, C):
    """
    PURPOSE
    -------
        Function calculates the change in energy particle due to collisions 
        in thermal plasma.

    INPUTS
    ------


    OUTPUTS
    -------
        de_dn: float
            Change in energy for that distance traveled.
    """
    de_dn = (-1*C)/(mu_naught*X*(1 - (3*C*n)/(mu_naught*X**2))**float(2/3)) 

    return de_dn

def calculate_change_in_energy(energy_value, dn, electron_nd=1e9*u.cm**-3, ion_nd=1e9*u.cm**-3, Te=2e6*u.K, Ti=2e6*u.K):
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

        electron_nd : float (optional)
            Electron number density

        ion_nd : float (optional)
            Ion number density

        Te : float (optional)
            Electron temperature

        Ti : float (optional)
            Ion temperature 

    OUTPUTS
    -------
        n: list
            List of the energies of the particle at each distance step
    """
    energy_J = energy_value.to((u.kg*u.m**2)/u.s**2)
    velocity_mag = np.sqrt(2*energy_J/E_MASS_G.to(u.kg))
    vx, vy, vz = np.sqrt((velocity_mag**2)/3), np.sqrt((velocity_mag**2)/3), np.sqrt((velocity_mag**2)/3)
    energy_value = energy_value.to(u.erg)
    mu_naught = vz/velocity_mag
    coulomb_log = calculate_coulomb_log(electron_nd=electron_nd, ion_nd=ion_nd, Te=Te, Ti=Ti)
    C = 2*np.pi*(E_CHARGE_ESMU**4)*coulomb_log

    stop_distance = collisions(energy_value, [vx, vy, vz])
    ns = np.arange(0, stop_distance.value, dn)
    new_energy = integrate.odeint(change_in_energy, energy_value.value, ns, args=(mu_naught, C))
    new_energy = new_energy * u.erg
    new_energy = new_energy.to(u.keV)
    e_thermal = float(3/2)*k_B*Te
    e_thermal = e_thermal.to(u.keV)

    return new_energy
