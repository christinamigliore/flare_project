"""
PURPOSE
-------
	Calculates the position and velocity x, y, z components of
    charged particles moving through a B field and E field using 
    the Lorentizan force and the kinematic equations of motion.
    Additionally 3D plots the trajectories of the particle. 
    IN SI UNITS!
"""
import numpy as np
import unittest
import multiprocessing
from scipy import integrate

import random_generator as rg
import stats
import velocity_particles as vp

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

# CONSTANTS
# Electron mass in kilograms
E_MASS_G = 9.1094e-31
# Proton mass in kilograms
P_MASS_G = 1.6726e-27
# Speed of light [m/s]
SPEED_OF_LIGHT = 2.9979e8
ELECTRON_CHARGE = 1.6021766208e-19


def solver(X, t0, Efield, Bfield, mass_of_particle):
    """
    PURPOSE:
    --------
        Function takes in the mass of particle, the initial velocity, the initial position, 
        the B field, and E field from the motion_parallel or motion_not_parallel function. It defines
        new variables for each x, y, z components of these inputs and calculates the lorentzian 
        acceleration. This function is then integrated using ODE solver in scipy.

    INPUTS
    ------
        X: Nd array
            X[1] = x (intial x postion)
            X[2] = y (intial y postion)
            X[3] = z (intial z postion)
            X[4] = Vx (intial x velocity)
            X[5] = Vy (intial y velocity)
            X[6] = Vz (intial z velocity)

        t0 = 1d array
            Array of time points

        B_field: array
            The magnetic field. [Bx, By, Bz]

        E_field: array
            The electric field. [Ex, Ey, Ez]

        mass_of_particle: float
            Mass of the particle

    OUTPUTS
    -------
        [vx, vy, vz, ax, ay, az]: array
            Array of calculated acceleration and velocity of particle, by components. Used for the integration
            over all the time points.
    """
    x = X[0]
    y = X[1]
    z = X[2]

    vx = X[3]
    vy = X[4]
    vz = X[5]

    ax = ELECTRON_CHARGE*(Efield[0]+(vy*Bfield[2])-(vz*Bfield[1]))/mass_of_particle
    ay = ELECTRON_CHARGE*(Efield[1]+(vz*Bfield[0])-(vx*Bfield[2]))/mass_of_particle
    az = ELECTRON_CHARGE*(Efield[2]+(vx*Bfield[1])-(vy*Bfield[0]))/mass_of_particle
    return [vx, vy, vz, ax, ay, az ]

def motion_not_parallel(v_init, p_init, B_field=[0, 0, 0.1], E_field=[0, 0, 0], ticks=2100, t_length=1e-8, particle_type='electron'):
    """
    PURPOSE:
    --------
        The function assigns a mass depending on the inputed particle 
        then sets up a grid space using the particle's gyrodiameter as the domain. The
        grid is then divided into the specified number of ticks, from this the time step
        is calculated. Using the time array created, which is set to sample twice at each
        tick in the grid, the function calls upon the solver function to find the particle's
        velocity and position at each sampling point. 

    INPUTS:
    ------
        v_init: array
            The initial velocity of particle in component form. [Vx0, Vy0, Vz0]

        p_init: array
            The initial positon of particle. [x0, y0, z0]

        B_field: array (optional)
            The magnetic field. [Bx, By, Bz]

        E_field: array (optional)
            The electric field. [Ex, Ey, Ez]

        ticks: integer (optional)
            Number of times for divide domain into (domain = gyrodiameter of particle)

        t_length: float (optional)
            The duration the particle travels in this simulated grid

        particle_type: string (optional)
            Can be either proton or electron. Set to electron as default.
            
    OUTPUTS:
    -------
        result: Nd array
            Nd array list of particle's positions and velocities
    """
    if particle_type == 'proton':
        mass_of_particle = P_MASS_G
    else:
        mass_of_particle = E_MASS_G

    velocity_mag = np.sqrt(v_init[0]**2 + v_init[1]**2 + v_init[2]**2)
    gyro_diameter = (2*(mass_of_particle*velocity_mag)/(ELECTRON_CHARGE*B_field[2]))
    tick_length = gyro_diameter/ticks
    time_step = tick_length/velocity_mag
    print time_step

    t = np.arange(0, t_length, time_step/2)
    pv0 = [p_init[0], p_init[1], p_init[2], v_init[0], v_init[1], v_init[2]]

    pv = integrate.odeint(solver, pv0, t, args=(E_field, B_field, mass_of_particle))
    return pv 


def motion_parallel(params):
    """
    PURPOSE:
    --------
        The function assigns a mass depending on the inputed particle 
        then sets up a grid space using the particle's gyrodiameter as the domain. The
        grid is then divided into the specified number of ticks, from this the time step
        is calculated. Using the time array created, which is set to sample twice at each
        tick in the grid, the function calls upon the solver function to find the particle's
        velocity and position at each sampling point. This function is parallelized and therefore
        can take more than one particle at a time, run in conjuction with the function "run_particles"
        which specifies the number of particle to simulate.

    INPUTS:
    ------
        params: Nd array
            params[0] = v_init, params[1] = p_init, params[2] = B_field, parmas[3] = E_field, params[4] = grid 
            where grid[0] = ticks, grid[1] = t_length, grid[2] = particle_type

            v_init: array
                The initial velocity of particle in component form. [Vx0, Vy0, Vz0]

            p_init: array
                The initial positon of particle. [x0, y0, z0]

            B_field: array (optional)
                The magnetic field. [Bx, By, Bz]

            E_field: array (optional)
                The electric field. [Ex, Ey, Ez]

            ticks: integer (optional)
                Number of times for divide domain into (domain = gyrodiameter of particle)

            t_length: float (optional)
                The duration the particle travels in this simulated grid

            particle_type: string (optional)
                Can be either proton or electron. Set to electron as default.
            
    OUTPUTS:
    -------
         result: Nd array
            Nd array list of particle's positions and velocities 
    """
    v_init, p_init, B_field, E_field, grid = params
    ticks, t_length, particle_type = int(grid[0]), float(grid[1]), grid[2]

    if particle_type == 'proton':
        mass_of_particle = P_MASS_G
    else:
        mass_of_particle = E_MASS_G

    velocity_mag = np.sqrt(v_init[0]**2 + v_init[1]**2 + v_init[2]**2)
    gyro_diameter = (2*(mass_of_particle*velocity_mag)/(ELECTRON_CHARGE*B_field[2]))
    tick_length = gyro_diameter/ticks
    time_step = tick_length/velocity_mag

    t = np.arange(0, t_length, time_step/2)
    pv0 = [p_init[0], p_init[1], p_init[2], v_init[0], v_init[1], v_init[2]]

    pv = integrate.odeint(solver, pv0, t, args=(E_field, B_field, mass_of_particle))
    return pv 

def run_particles(num_samples, t_length, Bfield=[0, 0, 0.1], Efield=[0, 0, 0], ticks=2100, particle_type='electron'):
    """
    PURPOSE:
    --------
        The function assigns a mass depending on the inputed particle 
        then sets up a grid space using the particle's gyrodiameter as the domain. The
        grid is then divided into the specified number of ticks, from this the time step
        is calculated. Using the time array created, which is set to sample twice at each
        tick in the grid, the function calls upon the solver function to find the particle's
        velocity and position at each sampling point. This function is parallelized and therefore
        can take more than one particle at a time, run in conjuction with the function "run_particles"
        which specifies the number of particle to simulate.

    INPUTS:
    ------
        num_samples: integer
            Number of particles to simulate

        t_length: float (optional)
            The duration the particle travels in this simulated grid

        Bfield: array (optional)
            The magnetic field. [Bx, By, Bz]

        Efield: array (optional)
            The electric field. [Ex, Ey, Ez]

        ticks: integer (optional)
            Number of times for divide domain into (domain = gyrodiameter of particle)

        particle_type: string (optional)
            Can be either proton or electron. Set to electron as default.
            
    OUTPUTS:
    -------
         result: Nd array
            Nd array list of particle's positions and velocities 
    """
    x_vals = stats.sample_cdf(num_samples)
    #vel = vp.get_vx_vy_vz(x_vals, len(x_vals))
    vel = rg.get_random_points_velocity(len(x_vals), 1, 3, 3)
    pos = rg.get_random_points(len(x_vals), 0, 1, 3)
    vel[0] = [3e5, 3e5, 3e5]
    pos[0] = [0, 0, 0]
    Bfields = np.full((len(vel), 3), Bfield)
    Efields = np.full((len(vel), 3), Efield)
    grid_params = np.full((len(vel), 3), [ticks, t_length, particle_type])
    params = zip(vel, pos, Bfields, Efields, grid_params)
    p = multiprocessing.Pool()
    result = p.map(motion_parallel, params)

    #plot(result)

    return result

def plot(result):
    """
    PURPOSE:
    --------
        This function runs finding_velocity_and_position and returns x, y, z
        component position arrays for the number of particles specified and 
        3D plots them in parallel.

    INPUTS:
    -------
        result: Nd array
            result[particle number][time array, x y z vx vy or vz]
            index x = 0, y = 1, z = 2, vx = 3, vy = 4, vz = 5.

    OUTPUTS:
    --------
        Plot of the trajectories for each simulated particle.
    """
 
    mpl.rcParams['legend.fontsize'] = 10
    ax = Axes3D(figure())
    for value in result:
        ax.plot(value[:,0], value[:,1], value[:,2])  
    ax.set_zlabel('Z')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Particle Trajectory')
    show()

def test_ticks(ticks=100):
    """
    PURPOSE:
    --------
        This function tests the percent difference of the calculated velocities of the particle
        at different number of ticks at an interval of 100 ticks. Once the percent difference is
        less than 1, the code prints the number of ticks where this occurred.

    INPUTS:
    -------
        ticks: integer (optional)
            Is the intial amount of ticks to test at

    OUTPUTS:
    --------
        Prints the tick number that produced a percent difference of less than 1%
    """
    next_tick = ticks + 100
    pv1 = motion_not_parallel([1e5, 1e5, 1e5], [0, 0, 0], ticks=ticks)
    pv2 = motion_not_parallel([1e5, 1e5, 1e5], [0, 0, 0], ticks=next_tick)
    test_val_vel = [np.abs(((pv1[10, 3] - pv2[10, 3])/pv2[10, 3]))*100, np.abs(((pv1[10, 4] - pv2[10, 3])/pv2[10, 4]))*100, np.abs(((pv1[10, 5] - pv2[10, 5])/pv2[10, 5]))*100]
    test_val_pos = [np.abs(((pv1[10, 0] - pv2[10, 0])/pv2[10, 0]))*100, np.abs(((pv1[10, 1] - pv2[10, 1])/pv2[10, 1]))*100, np.abs(((pv1[10, 2] - pv2[10, 2])/pv2[10, 2]))*100]
    if test_val_vel[0] <= 1 and test_val_vel[1] <= 1 and test_val_vel[2] <= 1:
        print "NEW TICK VALUE: "+str(next_tick)
    else:
        print next_tick
        test_ticks(next_tick)

def test_energy(result):
    """
    PURPOSE:
    --------
        This function tests the energy of the particle as it moves through the grid.
        Since magnetic fields do not work, there should be energy conservation, 
        (1/2)mVi^2 = (1/2)mVf^2. This code tests this.

    INPUTS:
    -------
        result: Nd array
            Output of the motion_not_parallel or motion_parallel function (particle's velocities
            and positions at each point)

    OUTPUTS:
    --------
        energy: 1d array
            Total energy of the particle at each sampled point in its tragectory in eV. 
    """
    for value in result:
        vel = np.sqrt(value[:, 3]**2 + value[:, 4]**2 + value[:, 5]**2)
        energy = ((0.5)*E_MASS_G*vel**2)*6.242e+18
        print "Percent Diff of Vi and Vf: "+str(np.abs(energy[0] - energy[-1])*100/energy[0])
    return energy
