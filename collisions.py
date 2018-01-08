"""
Collisions module
IN CGS UNITS!
"""
import math
import numpy as np
from astropy import units as u
from scipy import integrate
import matplotlib.pyplot as plt

import debye
import get_b_crit
import energy2vel as e
import column_depth as cd

# Charge of electron in esmu
E_CHARGE = 4.80e-10
# Electron mass in grams
E_MASS_G = 9.1094e-28*u.g
# Proton mass in grams
P_MASS_G = 1.6726219e-24*u.g

def change_in_energy(t, y):
    """
    PURPOSE 
    -------
        This function calculates the change in energy of a particle from collisions.
        This is an inside function meant for an integration solver.

    INPUTS
    ------
        t : float
            time 

        y : float
            energy at time t

    OUTPUTS
    -------
        dE_dt : float 
            Instantaneous change in energy at time t
    """
    coulomb_log = 20.825
    density = 1e9

    energy = y[0]
    n = len(y) 
    dE_dt = np.zeros((n,1))
    dE_dt[0] = ((-2*np.pi*(E_CHARGE**4)*coulomb_log)/(energy))*(E_MASS_G/P_MASS_G)*density*np.sqrt((2*energy)/E_MASS_G)

    return dE_dt

def change_in_vel(t, y):
    """
    PURPOSE 
    -------
        This function calculates the change in velocity (parallel) of a particle from collisions.
        This is an inside function meant for an integration solver.

    INPUTS
    ------
        t : float
            time 

        y : float
            energy at time t

    OUTPUTS
    -------
        dV_dt : float 
            Instantaneous change in velocity (parallel) at time t
    """
    coulomb_log = 20.825
    density = 1e9
    vz = 10828409511.3

    velocity_parallel = y[0]

    n = len(y) 
    dV_dt = np.zeros((n,1))
    dV_dt[0] = ((-1*np.pi*(E_CHARGE**4)*coulomb_log)/((0.5*E_MASS_G*(velocity_parallel**2 + vz**2))*2))*(1 + (E_MASS_G/P_MASS_G))*(velocity_parallel**2 + vz**2)

    return dV_dt


def collisions(energy, ne=1e9*u.cm**-3, ni=1e9*u.cm**-3, Te=2e6*u.K, Ti=2e6*u.K):
    """
    PURPOSE
    -------
        This function calculates the change in energy and parallel velocity from particle
        collisions. 

    INPUTS
    ------
        energy : float * astropy unit
            Initial energy of particle (can be in any units)

        ne : float (optional)
            Electron number density

        ni : float (optional)
            Ion number density

        Te : float (optional)
            Electron temperature

        Ti : float (optional)
            Ion temperature

    OUTPUTS
    -------
        E : array
            The energy of the particle at each instant in time t

        V_par : array
            The parallel velocity of the particle at each instant in time t

    """
    energy_J = energy.to((u.kg*u.m**2)/u.s**2)
    velocity_mag = np.sqrt(2*energy_J/E_MASS_G.to(u.kg))
    vx, vy, vz = np.sqrt((velocity_mag**2)/3), np.sqrt((velocity_mag**2)/3), np.sqrt((velocity_mag**2)/3)
    vx = vx.to(u.cm/u.s)
    vy = vy.to(u.cm/u.s)
    vz = vz.to(u.cm/u.s)
    velocity_mag = velocity_mag.to(u.cm/u.s)
    e0_erg = energy_J.to(u.erg)
    velocity_parallel = np.sqrt(vy**2 + vz**2)
    print('PARAMETERS:')
    print("-----------------------------")
    print("Energy in keV: "+str(e0_erg.value*6.242e8))

    coulomb_log = cd.calculate_coulomb_log(electron_nd=ne, ion_nd=ni, Te=Te, Ti=Ti)
    print("Coulomb log: "+str(coulomb_log))

    t_start = 0.0
    t_end = 250000.0
    dt = 1e3
    num_steps = int(np.floor((t_end - t_start)/dt) + 1)

    r1 = integrate.ode(change_in_energy).set_integrator('vode', method='bdf')
    r1.set_initial_value([e0_erg.value], t_start)
    t = np.zeros((num_steps, 1))
    E = np.zeros((num_steps, 1))
    t[0] = t_start
    E[0] = e0_erg.value * 6.242e8

    k = 1
    while r1.successful() and k < num_steps:
        r1.integrate(r1.t + dt)
        t[k] = r1.t
        if r1.y[0] * 6.242e8 > 0:
            E[k] = r1.y[0] * 6.242e8
        else:
            break 
        k += 1

    r2 = integrate.ode(change_in_vel).set_integrator('vode', method='bdf')
    r2.set_initial_value([velocity_parallel.value], t_start)
    V_par = np.zeros((num_steps, 1))
    V_par[0] = velocity_parallel.value

    j = 1
    while r2.successful() and j < num_steps:
        r2.integrate(r2.t + dt)
        V_par[k] = r2.y[0]
        j += 1

    #E_flattened = np.ndarray.flatten(E)
    #idx = np.argwhere(E_flattened == 0)[0][0]
    plt.figure()
    plt.subplot(211)
    plt.plot(t, E)
    plt.xlabel('Time (s)')
    plt.ylabel('Energy (keV)')
    plt.subplot(212)
    plt.plot(t, V_par)
    plt.xlabel('Time (s)')
    plt.ylabel('Velocity (cm/s)')
    plt.tight_layout()
    plt.show()
    return E, V_par
