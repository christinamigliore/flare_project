"""
DOCUMENTATION

PURPOSE
-------
    Produces lists of randomly generated initial velocities and initial positions.

INPUTS
------
    num_of_particles: integer
        Number of particles to be generated
        
OUTPUTS
-------
    particle_positions: list
        List of initial x, y, and z velocities of particles.

    particle_velocities: list
        List of initial x, y, and z positions of particles.
"""
import numpy as np

def get_particle_points(num_of_particles):

    position_vector = []
    velocity_vector = []

    for i in range(0, num_of_particles):
        position = np.random.uniform(0, 1, 3)
        velocity = np.random.uniform(10, 200, 3)
        position_vector.append(position)
        velocity_vector.append(velocity)

    return position_vector, velocity_vector