"""
Collisions module
IN CGS UNITS!
"""
import math
import numpy as np
from astropy import units as u
from scipy import integrate

import debye
import get_b_crit
import energy2vel as e

# Charge of electron in esmu
E_CHARGE = 4.80e-10
# Electron mass in grams
E_MASS_G = 9.1094e-28
# Proton mass in grams
P_MASS_G = 1.6726219e-24
# Electron density in cm^3
E_DENSITY = 1e9

def change_in_vel(X, t0, mu, debye_length):
    """
    """
    xi = 0.5*((E_CHARGE**2)/(mu*X**2))**2*np.log(1 + ((mu*debye_length.value*X**2)/E_CHARGE**2)**2)
    dv_dt = -1*4*np.pi*E_DENSITY*(P_MASS_G/(P_MASS_G+E_MASS_G))*xi*X**2

    return dv_dt

def change_in_energy(X, t0, mu, debye_length):
    """
    """
    xi = 0.5*((E_CHARGE**2)/(mu*X**2))**2*np.log(1 + ((mu*debye_length.value*X**2)/E_CHARGE**2)**2)
    dE_dt = 4*np.pi*E_DENSITY*((P_MASS_G*E_MASS_G**2)/(P_MASS_G+E_MASS_G)**2)*xi*X**3

    return dE_dt

def collisions(velocity, t_end=1, dt=0.01):
    """
    """

    vx = velocity[0] * u.m/u.s
    vy = velocity[1] * u.m/u.s
    vz = velocity[2] * u.m/u.s

    vx = vx.to(u.cm/u.s)
    vy = vy.to(u.cm/u.s)
    vz = vz.to(u.cm/u.s)

    velocity_mag = np.sqrt(vx.value**2 + vy.value**2 + vz.value**2)

    debye_length = debye.calculate_debye()
    debye_length = debye_length.to(u.cm)

    mu = (E_MASS_G*P_MASS_G)/(E_MASS_G+P_MASS_G)
    t0 = 0.1
    t = np.arange(t0, t_end, dt)
    v = integrate.odeint(change_in_vel, velocity_mag, t, args=(mu, debye_length), mxstep=5000)
    
    #r = integrate.ode(change_in_vel).set_integrator('vode', nsteps=5000, method='bdf')
    #r.set_initial_value(velocity_mag, t0).set_f_params(mu, debye_length)
    #while r.successful() and r.t < t_end:
        #print(r.t+dt, r.integrate(r.t+dt))

    return v

