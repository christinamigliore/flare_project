"""
PURPOSE
-------
    Time testing of the particle_motion_b_field. Calculates time between
    particle simulations then plots as a function of number of particles
    and runtime.
"""
import datetime as d
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np

import stats
import random_generator as rg 
import particle_motion_b_field as pb 
import particle_parallel as pp
import velocity_particles as vp

def calculate_delta_t(particle_num, B_field):
    """
    PURPOSE:
    --------
        This functions runs the module which calculates the positions and velocities 
        for the number of specified simulated particles and times it.

    INPUTS:
    -------
        particle_num: integer
            The number of simulated particles to produce.

        B_field: array
            [Bx, By, Bz], the magnetic field used.

    OUTPUTS:
    --------
        delta_time: float
            The time it takes to simulate the specified particle's velocitites and 
            positions.
    """
    time_start = d.datetime.now()
    x_vals = stats.sample_cdf(particle_num)
    vel = vp.get_vx_vy_vz(x_vals, len(x_vals))
    pos = rg.get_random_points(len(x_vals), 0, 1, 3)

    for i in range(0, len(vel)):
        v_init = vel[i] 
        p_init = pos[i] 
        particle_velocity, particle_position, length_of_grid = pb.equations_of_motion(v_init, p_init, B_field, particle_type='electron')
    time_end = d.datetime.now()
    delta_time = (time_end - time_start).total_seconds()

    print str(particle_num)+" particles, loop time: "+str(delta_time)
    return delta_time

def calculate_delta_t_parallel(particle_num, B_field):
    """
    PURPOSE:
    --------
        This functions runs the module which calculates the positions and velocities 
        for the number of specified simulated particles in parallel and times it.

    INPUTS:
    -------
        particle_num: integer
            The number of simulated particles to produce.

        B_field: array
            [Bx, By, Bz], the magnetic field used.

    OUTPUTS:
    --------
        delta_time: float
            The time it takes to simulate the specified particle's velocitites and 
            positions.
    """
    time_start = d.datetime.now()
    x_vals = stats.sample_cdf(particle_num)
    vel = vp.get_vx_vy_vz(x_vals, len(x_vals))
    pos = rg.get_random_points(len(x_vals), 0, 1, 3)
    B_fields = np.full((len(vel), 3), B_field)
    params = zip(vel, pos, B_fields)
    p = multiprocessing.Pool()
    p.map(pp.equations_of_motion, params)

    time_end = d.datetime.now()
    delta_time_p = (time_end - time_start).total_seconds()

    print str(particle_num)+" particles, parallel, loop time: "+str(delta_time_p)
    return delta_time_p

if __name__ == '__main__':
    num_of_particles_string = raw_input('Type particles numbers, ex. 1 10 100 ...: ')
    num_of_particles_list = num_of_particles_string.split()
    num_of_particles = [int(a) for a in num_of_particles_list]
    B_field_string = raw_input('Magnetic field (component-wise), ex. 0 0 1: ')
    B_field_list = B_field_string.split()
    B_field = [int(a) for a in B_field_list]
    delta_times = []
    for i in range(len(num_of_particles)):
    	delta_time = calculate_delta_t(num_of_particles[i], B_field)
    	delta_times.append(delta_time)

    delta_times_parallel = []
    for i in range(0, len(num_of_particles)):
    	delta_time_p = calculate_delta_t_parallel(num_of_particles[i], B_field)
    	delta_times_parallel.append(delta_time_p)

    plt.plot(num_of_particles, delta_times, label="non-parallel")
    plt.scatter(num_of_particles, delta_times, label="non-parallel points")
    plt.plot( num_of_particles, delta_times_parallel, label="parallel")
    plt.scatter( num_of_particles, delta_times_parallel, label="parallel points")
    plt.legend(loc=2)
    plt.ylabel("delta time")
    plt.xlabel("number of simulated particles")
    plt.ylim(0,)
    plt.xlim(0,)
    plt.show()
