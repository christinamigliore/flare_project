"""
PURPOSE
-------
    Calculates the position and velocity x, y, z components of
    charged particles moving through a uniform B field using 
    the Lorentizan force and the kinematic equations of motion.
    Additionally 3D plots the trajectories of each particle
    on the same plot. 
"""
import numpy as np
import unittest

import random_generator as rg
import stats
import velocity_particles as vp

# CONSTANTS
# Electron mass in grams
E_MASS_G = 9.1094*10**-28
# Proton mass in grams
P_MASS_G = 1.6726*10**-24
# Speed of light [cm/s]
SPEED_OF_LIGHT = 2.9979*10**10
ELECTRON_CHARGE = 1.6021766208 * 10**-19


def finding_velocity_and_position(mass_of_particle, v_next, next_position, B_field, time_step):
    """
    PURPOSE:
    ---------
	    Function takes in the mass of particle, the initial velocity, the initial position, 
	    the B field, and the calculated time step from the equations_of_motion function. It defines
	    new variables for each x, y, z components of these inputs and calculates the lorentzian 
	    acceleration. This acceleration is then plugged into the kinetic equations of motion
	    to obtain the final velocity and position. The x, y, z velocities and positions are then
	    each converted into an array. 

	INPUTS:
	-------
		mass_of_particle: float
			Mass of particle either electron or proton.

		v_next: array
			Array of the x, y, z velocity components of particle at one instance in time.

		next_position: array
			Array of the x, y, z position components of particle at one instance in time.

		B_field: array
			Array of the x, y, z magnetic field.

		time_step: array
			The time step for the area of the box. [start, end]

	OUTPUTS:
	--------
		velocity_vector: array
			Particle's next x, y, z velocity components

		position_vector: array
			Particle's next x, y, z position components
    """
    v_x = v_next[0]
    v_y = v_next[1]
    v_z = v_next[2]
    
    Bx = B_field[0]
    By = B_field[1]
    Bz = B_field[2]

    p_x = next_position[0]
    p_y = next_position[1]
    p_z = next_position[2]

    ax = (ELECTRON_CHARGE/mass_of_particle)*(v_y*Bz - v_z*By)/SPEED_OF_LIGHT
    ay = (ELECTRON_CHARGE/mass_of_particle)*(v_z*Bx - v_x*Bz)/SPEED_OF_LIGHT
    az = (ELECTRON_CHARGE/mass_of_particle)*(v_x*By - v_y*Bx)/SPEED_OF_LIGHT


    time_step = time_step[1] - time_step[0]

    vx_final = v_x + ax*time_step
    vy_final = v_y + ay*time_step
    vz_final = v_z + az*time_step

    x_final = p_x + v_x*time_step + (1/2)*ax*time_step**2
    y_final = p_y + v_y*time_step + (1/2)*ay*time_step**2
    z_final = p_z + v_z*time_step + (1/2)*az*time_step**2

    position_vector = np.array([x_final, y_final, z_final])
    velocity_vector = np.array([vx_final, vy_final, vz_final])
    position_vector = position_vector.flatten()
    velocity_vector = velocity_vector.flatten()
    return velocity_vector, position_vector

def equations_of_motion(v_init, p_init, B_field, particle_type='electron'):
    """
    PURPOSE:
    --------
	    This function takes in the initial velocity, the B field, and the length of the 
	    grid as inputs. The function assigns a mass depending on the inputed particle 
	    then creates an empty array that is twice as long as the length of the grid
	    for the velocity and position array. The function then calculates the time step
	    by taking the distance traveled and dividing by twice the initial velocity.
	    Using a loop, the code interates through 2*length of grid and defines a new time
	    step along with running the finding_velocity_and_position function for each 
	    point in the grid. The function then defines x_array, y_array, and z_array
	    in order to separate each of the position components in order to plot them.

	INPUTS
	------
		v_init_list: array
			The initial list of velocities for the particles. [[Vx10, Vy10, Vz10], ...]

		B_field: array
			The uniform magnetic field. [Bx, By, Bz]

		particle_type: string (optional)
			Can be either proton or electron. Set to electron as default.
			

	OUTPUTS
	-------
	    particle_velocity: 2D array
	        Contains all the calculated final velocities, Vx, Vy and Vz, of the 
	        particle at each point.

	    particle_position: 2D array
	        Containes all the calculated final positions, x, y and z, of the
	        particle at each point.

	    length_of_grid: integer
	        The width, height, length of the 3D grid box.
	"""
    if particle_type == 'proton':
        mass_of_particle = P_MASS_G
    else:
        mass_of_particle = E_MASS_G
        
    time_step = 1/float(2*(v_init[0]))
    characteristic_length = time_step*v_init[0]
    length_of_grid = int(5000*characteristic_length)  
    
    particle_velocity = np.zeros((2*length_of_grid, 3))
    particle_velocity[0] = v_init
    particle_position = np.zeros((2*length_of_grid, 3))
    particle_position[0] = p_init
    
    for i in range(0, 2*length_of_grid-1):
        new_time_step = [0+i, time_step+i]
        particle_velocity[i+1], particle_position[i+1] = finding_velocity_and_position(mass_of_particle, particle_velocity[i], particle_position[i], B_field, new_time_step)

    return particle_velocity, particle_position, length_of_grid

def run_particles(num_samples, B_field):
	"""
	PURPOSE:
	--------
		This functions runs the position and velocity analysis and
		returns the positions and velocities for the user inputted 
		number of particles.

	INPUTS:
	-------
		num_samples: integer
			Number of particles to simulate

		B_field: array
            The uniform magnetic field. [Bx, By, Bz]

	OUTPUTS:
	--------
		particle_velocities: list
			List of all the positions for the simluated particles

		particle_positions: list
			List of all the velocities for the simulated particles 
	"""
	x_vals = stats.sample_cdf(particle_num)
    vel = vp.get_vx_vy_vz(x_vals, len(x_vals))
    pos = rg.get_random_points(len(x_vals), 0, 1, 3)
	particle_velocities = []
	particle_positions = []
	for i in range(0, len(vel)):
		v_init = vel[i] 
		p_init = pos[i] 
		particle_velocity, particle_position, length_of_grid = equations_of_motion(v_init, p_init, B_field, particle_type='electron')
		particle_velocities.append(particle_velocity)
		particle_positions.append(particle_position)
	return particle_velocities, particle_positions
"""
class Test_finding_velocity_and_position(unittest.TestCase):
    def test_finding_velocity_and_position(self):
        self.assertEqual(finding_velocity_and_position(9.1094*10**-28, [100, 10, 30], [0, 0, 0], [0, 0, 10], [0, 0.5]), [100.293341481, 7.06658519396, 30], [50.2200061105, 2.79993889547, 15.0])

if __name__ == '__main__':
    unittest.main()
"""
