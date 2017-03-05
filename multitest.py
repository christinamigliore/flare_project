"""
DOCUMENTATION

PURPOSE
-------
    Calculates the position and velocity x, y, z components of
    charged particles moving through a uniform B field using 
    the Lorentizan force and the kinematic equations of motion.
    Additionally 3D plots the trajectories of each particle
    on the same plot. This is done through multiprocessing.

INPUTS
------
	v_init_list: array
		The initial list of velocities for the particles. [[Vx10, Vy10, Vz10], ...]

	B_field: array
		The uniform magnetic field. [Bx, By, Bz]

	particle_type: string
		Can be either proton or electron. Set to electron as default.
		

OUTPUTS
-------
    particle_velocity: 2D array
        Contains all the calculated final velocities, Vx, Vy and Vz, of the 
        particle at each point.

    particle_position: 2D array
        Containes all the calculated final positions, x, y and z, of the
        particle at each point.

    x_array: 1D array
        Array of all the x positions of the particle. 

    y_array: 1D array
        Array of all the y positions of the particle.

    z_array: 1D array
        Array of all the z positions of the particle. 

    length_of_grid: integer
        The width, height, length of the 3D grid box.

"""
import numpy as np
from scipy.integrate import odeint
import unittest
import multiprocessing 

import random_generator as rg

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
    Function takes in the mass of particle, the initial velocity, the initial position, 
    the B field, and the calculated time step from the equations_of_motion function. It defines
    new variables for each x, y, z components of these inputs and calculates the lorentzian 
    acceleration. This acceleration is then plugged into the kinetic equations of motion
    to obtain the final velocity and position. The x, y, z velocities and positions are then
    each converted into an array. 
    """
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
    
    def dvx_dt(vx, t):
        vx = (ELECTRON_CHARGE/mass_of_particle)*(v_y*Bz - v_z*By)/SPEED_OF_LIGHT
        return vx    
    def dvy_dt(vy, t):
        vy = (ELECTRON_CHARGE/mass_of_particle)*(v_z*Bx - v_x*Bz)/SPEED_OF_LIGHT
        return vy  
    def dvz_dt(vz, t):
        vz = (ELECTRON_CHARGE/mass_of_particle)*(v_x*By - v_y*Bx)/SPEED_OF_LIGHT
        return vz
    
    velocity_x = odeint(dvx_dt, v_x, time_step)
    velocity_y = odeint(dvy_dt, v_y, time_step)
    velocity_z = odeint(dvz_dt, v_z, time_step)
    
    def dx_dt(x, t):
        x = velocity_x[1]
        return x    
    def dy_dt(y, t):
        y = velocity_y[1]
        return y  
    def dz_dt(z, t):
        z = velocity_z[1]
        return z

    position_x = odeint(dx_dt, p_x, time_step)
    position_y = odeint(dy_dt, p_y, time_step)
    position_z = odeint(dz_dt, p_z, time_step)
    
    position_vector = np.array([position_x[1], position_y[1], position_z[1]])
    velocity_vector = np.array([velocity_x[1], velocity_y[1], velocity_z[1]])
    position_vector = position_vector.flatten()
    velocity_vector = velocity_vector.flatten()
    return velocity_vector, position_vector
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

def equations_of_motion(params, B_field=[0, 0, 100], particle_type='electron'):
    """
    This function takes in the initial velocity, the B field, and the length of the 
    grid as inputs. The function assigns a mass depending on the inputed particle 
    then creates an empty array that is twice as long as the length of the grid
    for the velocity and position array. The function then calculates the time step
    by taking the distance traveled and dividing by twice the initial velocity.
    Using a loop, the code interates through 2*length of grid and defines a new time
    step along with running the finding_velocity_and_position function for each 
    point in the grid. The function then defines x_array, y_array, and z_array
    in order to separate each of the position components in order to plot them.
    """
    v_init, p_init = params
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

if __name__ == '__main__':
    pos, vel = rg.get_particle_points(5)
    params = zip(vel, pos)
    p = multiprocessing.Pool()
    p.map(equations_of_motion, params)

"""
class Test_finding_velocity_and_position(unittest.TestCase):
    def test_finding_velocity_and_position(self):
        self.assertEqual(finding_velocity_and_position(9.1094*10**-28, [100, 10, 30], [0, 0, 0], [0, 0, 10], [0, 0.5]), [100.293341481, 7.06658519396, 30], [50.2200061105, 2.79993889547, 15.0])

if __name__ == '__main__':
    unittest.main()
"""
