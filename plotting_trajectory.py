"""
PURPOSE
-------
 	Function calculates the randomly generated particles'
 	trajectories and plots them in a 3D plot either using
    series calculations or parallel.
"""
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import matplotlib.animation as animation

import particle_motion_b_field as pm 
import particle_parallel as pp 
import multitest as mt
import random_generator as rg
import stats
import velocity_particles as vp

def plot_particle_motion(num_samples, B_field):
    """
    PURPOSE:
    --------
        This function runs finding_velocity_and_position and returns x, y, z
        component position arrays for each inputed initial velocity (i.e. particle)
        and 3D plots them. The axis are set to the length of grid.  

    INPUTS
    ------
        num_samples: integer
            Number of particles to simulate

        B_field: array
            The uniform magnetic field. [Bx, By, Bz]
    """
    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    x_vals = stats.sample_cdf(num_samples)
    v_init_list = vp.get_vx_vy_vz(x_vals, len(x_vals))
    p_init_list = rg.get_random_points(len(x_vals), 0, 1, 3)

    for i in range(0, len(v_init_list)):
        v_init = v_init_list[i] 
        p_init = p_init_list[i] 
        particle_velocity, particle_position, length_of_grid = pm.equations_of_motion(v_init, p_init, B_field, particle_type='electron')

        x_array = np.zeros(len(particle_position))
        y_array = np.zeros(len(particle_position))
        z_array = np.zeros(len(particle_position))
    
        for j in range(0, len(particle_position)):
            x_array[j] = particle_position[j][0]
            y_array[j] = particle_position[j][1]
            z_array[j] = particle_position[j][2]

        ax.plot(x_array, y_array, z_array)
        ax.set_xlim(0, length_of_grid)
        ax.set_ylim(0, length_of_grid)
        ax.set_zlim(0, length_of_grid)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
    plt.show()

def plot_particle_motion_parallel(num_samples, B_field):
    """
    PURPOSE:
    --------
        This function runs finding_velocity_and_position and returns x, y, z
        component position arrays for the number of particles specified and 
        3D plots them in parallel. The axis are set to the length of grid.  

    INPUTS:
    -------
        num_samples: integer
            Number of particles to simulate

        B_field: array
            The uniform magnetic field. [Bx, By, Bz]
    OUTPUTS:
    --------
        Plot of the trajectories for each simulated particle.
    """

    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
 
    x_vals = stats.sample_cdf(num_samples)
    vel = vp.get_vx_vy_vz(x_vals, len(x_vals))
    pos = rg.get_random_points(len(x_vals), 0, 1, 3)
    B_fields = np.full((len(vel), 3), B_field)
    params = zip(vel, pos, B_fields)
    p = multiprocessing.Pool()
    result = p.map(pp.equations_of_motion, params)
    
    for value in result:
        x_array = np.zeros(len(value[0]))
        y_array = np.zeros(len(value[0]))
        z_array = np.zeros(len(value[0]))

        length_of_grid = value[2]

        for j in range(0, len(value[1])):
            x_array[j] = value[1][j][0]
            y_array[j] = value[1][j][1]
            z_array[j] = value[1][j][2]

        ax.plot(x_array, y_array, z_array)
        ax.set_xlim(0, length_of_grid)
        ax.set_ylim(0, length_of_grid)
        ax.set_zlim(0, length_of_grid)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
    plt.show()