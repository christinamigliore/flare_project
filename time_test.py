"""
DOCUMENTATION

PURPOSE
-------
    Time testing of the particle_motion_b_field. 
    Calculates time between particle simulations 
    then plots as a function of number of particles
    and runtime.

INPUTS
------
    num_of_particles: list
    	List of number of particles to be simulated
    	and tested.
        
OUTPUTS
-------
    delta_time: list
        List of runtimes for particle simulations
"""
num_of_particles = [10, 100, 1000, 2000, 5000]

import datetime as d
import matplotlib.pyplot as plt
import multiprocessing

import random_generator as rg 
import particle_motion_b_field as pb 
import multitest as mt 

def calculate_delta_t(particle_num):
    time_start = d.datetime.now()
    pos, vel = rg.get_particle_points(particle_num)

    for i in range(0, len(vel)):
        v_init = vel[i] 
        p_init = pos[i] 
        particle_velocity, particle_position, length_of_grid = pb.equations_of_motion(v_init, p_init, [0, 0, 100], particle_type='electron')
    time_end = d.datetime.now()
    delta_time = (time_end - time_start).total_seconds()

    print str(particle_num)+" particles, loop time: "+str(delta_time)
    return delta_time

def calculate_delta_t_parallel(particle_num):
    time_start = d.datetime.now()
    pos, vel = rg.get_particle_points(particle_num)
    params = zip(vel, pos)
    p = multiprocessing.Pool()
    p.map(mt.equations_of_motion, params)

    time_end = d.datetime.now()
    delta_time_p = (time_end - time_start).total_seconds()

    print str(particle_num)+" particles, parallel, loop time: "+str(delta_time_p)
    return delta_time_p

delta_times = []
for i in range(len(num_of_particles)):
	delta_time = calculate_delta_t(num_of_particles[i])
	delta_times.append(delta_time)

delta_times_parallel = []
for i in range(0, len(num_of_particles)):
	delta_time_p = calculate_delta_t_parallel(num_of_particles[i])
	delta_times_parallel.append(delta_time_p)

plt.plot(delta_times, num_of_particles, label="non-parallel")
plt.plot(delta_times_parallel, num_of_particles, label="parallel")
plt.legend()
plt.xlabel("delta time")
plt.ylabel("number of simulated particles")
plt.ylim(0,)
plt.xlim(0,)
plt.show()
