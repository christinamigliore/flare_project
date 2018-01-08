"""
Collisions module
(INTEGRATION EQUATIONS IN CGS UNITS!)
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
    vx, vy = 10828409511.3, 10828409511.3

    vz = y[0]

    n = len(y) 
    dV_dt = np.zeros((n,1))
    dV_dt[0] = ((-1*np.pi*(E_CHARGE**4)*coulomb_log)/((0.5*E_MASS_G*(vx**2 + vy**2+ vz**2))*2))*(1 + (E_MASS_G/P_MASS_G))*(vx**2 + vy**2 + vz**2)

    return dV_dt

def internal_velocity_function(vz, t_start, num_steps, dt):
    """
    PURPOSE 
    -------
        This function is an internal function to the function "collisions." This function uses integration
        solver to solve for the changes in parallel velocity.

    INPUTS
    ------
        vz : float
            The initial parallel velocity calculated in the function "collisions"

        t_start : float
            Start time

        num_steps : float
            Number of time num_steps

        dt : float
            The change of time in each step

    OUTPUTS
    -------
        V_par : array
            The parallel velocity of the particle at each instant in time t

        t : array
            The time array
    """
    r = integrate.ode(change_in_vel).set_integrator('vode', method='bdf')
    r.set_initial_value([vz.value], t_start)
    t = np.zeros((num_steps, 1))
    V_par = np.zeros((num_steps, 1))
    t[0] = t_start
    V_par[0] = vz.value

    k = 1
    while r.successful() and k < num_steps:
        r.integrate(r.t + dt)
        t[k] = r.t
        V_par[k] = r.y[0]
        k += 1
    return V_par, t


def internal_energy_function(e0_erg, t_start, num_steps, dt):
    """
    PURPOSE 
    -------
        This function is an internal function to the function "collisions." This function uses integration
        solver to solve for the changes in energy.

    INPUTS
    ------
        e0_erg : float
            The initial energy of the particle

        t_start : float
            Start time

        num_steps : float
            Number of time num_steps

        dt : float
            The change of time in each step

    OUTPUTS
    -------
        E : array
            The energy of the particle at each instant in time tx
    """
    r = integrate.ode(change_in_energy).set_integrator('vode', method='bdf')
    r.set_initial_value([e0_erg.value], t_start)
    t = np.zeros((num_steps, 1))
    E = np.zeros((num_steps, 1))
    t[0] = t_start
    E[0] = e0_erg.value * 6.242e8

    k = 1
    while r.successful() and k < num_steps:
        r.integrate(r.t + dt)
        t[k] = r.t
        if r.y[0] * 6.242e8 > 0:
            E[k] = r.y[0] * 6.242e8
        else:
            break 
        k += 1
    return E


def collisions(energy, vx, vy, vz, ne=1e9*u.cm**-3, ni=1e9*u.cm**-3, Te=2e6*u.K, Ti=2e6*u.K, t_start=0.0, t_end=250000.0, dt=1e3):
    """
    PURPOSE
    -------
        This function calculates the change in energy and parallel velocity from particle
        collisions. 

    INPUTS
    ------
        energy : float * astropy unit
            Initial energy of particle (can be in any units)
            Must include astropy unit!

        vx : float * astropy unit
            The x component of the particle's velocity (can be in any units)
            Must include astropy unit!

        vy : float * astropy unit
            The y component of the particle's velocity (can be in any units)
            Must include astropy unit!

        vz : float * astropy unit
            The z component of the particle's velocity (can be in any units)
            Must include astropy unit!

        ne : float (optional)
            Electron number density

        ni : float (optional)
            Ion number density

        Te : float (optional)
            Electron temperature

        Ti : float (optional)
            Ion temperature

        t_start : float (optional)
            Start time of simulation

        t_end : float (optional)
            End time of simulation

        dt : float (optional)
            The time step

    OUTPUTS
    -------
        E : array
            The energy of the particle at each instant in time t

        V_par : array
            The parallel velocity of the particle at each instant in time t

    """
    energy_J = energy.to((u.kg*u.m**2)/u.s**2)
    #velocity_mag = np.sqrt(2*energy_J/E_MASS_G.to(u.kg))
    #velocity_mag = velocity_mag.to(u.cm/u.s)
    #vx, vy, vz = np.sqrt((velocity_mag**2)/3), np.sqrt((velocity_mag**2)/3), np.sqrt((velocity_mag**2)/3)
    vx = vx.to(u.cm/u.s)
    vy = vy.to(u.cm/u.s)
    vz = vz.to(u.cm/u.s)
    e0_erg = energy_J.to(u.erg)
    print('PARAMETERS:')
    print("-----------------------------")
    print("Energy in keV: "+str(e0_erg.value*6.242e8))

    coulomb_log = cd.calculate_coulomb_log(electron_nd=ne, ion_nd=ni, Te=Te, Ti=Ti)
    print("Coulomb log: "+str(coulomb_log))

    num_steps = int(np.floor((t_end - t_start)/dt) + 1)

    E = internal_energy_function(e0_erg, t_start, num_steps, dt)
    V_par, t = internal_velocity_function(vz, t_start, num_steps, dt)

    plt.figure()
    plt.subplot(211)
    plt.title('Changes Due Collisions')
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
