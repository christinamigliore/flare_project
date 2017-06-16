import math
import numpy as np

import random_generator as rg 
import stats
import particle_motion_b_field as pb 
import plotting_trajectory as pt

def get_velocity_components():
    """
    PURPOSE:
    -------
        This function generates a random number from 0 to 1 for each x, y, z velocity
        components, representing the percentages of the total velocity in that direction
        and returns the total sum. If the total sum is zero then the function reruns to get
        a non zero value.

    INPUTS:
    ------
        None

    OUTPUTS:
    -------
        sum_vel: float
            The total sum of the x, y, z velocity components

        vx_precent: float
            The x velocity randomly generated value

        vy_precent: float
            The y velocity randomly generated value

        vz_present: float
            The z velocity randomly generated value
    """
    vx_percent = np.random.uniform(0, 1, 1)
    vy_percent = np.random.uniform(0, 1, 1)
    vz_percent = np.random.uniform(0, 1, 1)
    sum_vel = vx_percent + vy_percent + vz_percent
    sum_vel = float(sum_vel)
    if sum_vel != 0:
        return sum_vel, vx_percent, vy_percent, vz_percent
    else:
        sum_vel, vx_percent, vy_percent, vz_percent = get_velocity_components()
        return sum_vel, vx_percent, vy_percent, vz_percent

def percentage_vel(x_values, index):
    """
    PURPOSE:
    --------
        This function normalizes the randomly generated percentage velocity components
        found through the function get_velocity_components to have a sum of 1.
        The function then tests whether the new velocity components are a sum of 1
        then finds the actual velocities by using the position values. If the sum is not
        equal to 1, then the functions runs again.

    INPUTS:
    ------
        x_values: array
            Array of inputted position values for the particle.

        index: integer
            The index indicating which particle's velocity it is calculating.

    OUTPUTS:
    --------
        velocity_vector: array
            [Vx, Vy, Vz] initial position of simulated particle.
    """
    sum_vel, vx_percent, vy_percent, vz_percent = get_velocity_components()
    new_vx_percent = float(vx_percent)/sum_vel
    new_vy_percent = float(vy_percent)/sum_vel
    new_vz_percent = float(vz_percent)/sum_vel
    percent_sum = new_vx_percent + new_vy_percent + new_vz_percent
    if percent_sum == 1.0:
        vx = math.sqrt(x_values[index]*new_vx_percent)
        vy = math.sqrt(x_values[index]*new_vy_percent)
        vz = math.sqrt(x_values[index]*new_vz_percent)
        velocity_vector = np.array([vx, vy, vz])
        velocity_vector = velocity_vector.flatten()
        return velocity_vector
    else:
        velocity_vector = percentage_vel(x_values, index)
        return velocity_vector

def get_vx_vy_vz(x_values, num_samples):
    """
    PURPOSE:
    --------
        This functions randomly simulates the initial velocities for the user
        specified number of simulated particles using randomly simulated initial
        positions in the module particle_parallel.py or particle_motion_b_field.py.

    INPUTS:
    ------
        x_values: array
            Array of inputted position values for the particle.

        num_samples: integer
            Number of particles to simulate.

    OUTPUTS:
    --------
        velocity: array
            [[Vx1, Vy1, Vz1], ...] initial positions of all simulated particles.
    """
    velocity = np.zeros((num_samples, 3))
    for i in range(0, num_samples):
        velocity_vector = percentage_vel(x_values, i)
        velocity[i] = velocity_vector
    return velocity