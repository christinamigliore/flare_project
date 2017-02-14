"""
DOCUMENTATION

PURPOSE
-------
	Calculates the characteristic length scale of the inputed grid
	with the inputted parameter, (q).

INPUTS
------
	v_init: array
		The initial velocity of the particle. [Vx0, Vy0, Vz0]

	B_field: array
		The uniform magnetic field. [Bx, By, Bz]

	length_of_grid: integer
		The width, height, length of the 3D grid box.

	particle_type: string
		Can be either proton or electron. Set to electron as default.
		

OUTPUTS
-------
	particle_positons: array
		particle's postion at every time step.

"""
import numpy as np
from scipy.integrate import odeint
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

# CONSTANTS
# Electron mass in grams
E_MASS_G = 9.1094*10**-28
# Proton mass in grams
P_MASS_G = 1.6726*10**-24
# Speed of light [cm/s]
SPEED_OF_LIGHT = 2.9979*10**10
ELECTRON_CHARGE = 1.6021766208 * 10**-19

def finding_velocity_and_position(mass_of_particle, v_next, next_position, B_field, time_step):
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

def equations_of_motion(v_init, B_field, length_of_grid, particle_type='electron'):
    if particle_type == 'proton':
        mass_of_particle = P_MASS_G
    else:
        mass_of_particle = E_MASS_G

    particle_velocity = np.zeros((2*length_of_grid, 3))
    particle_velocity[0] = v_init
    particle_position = np.zeros((2*length_of_grid, 3))

    grid_dimensions = np.arange(length_of_grid)
    time_step = (grid_dimensions[1] - grid_dimensions[0])/float(2*(v_init[0]))

    for i in range(0, 2*length_of_grid-1):
        new_time_step = [0+i, time_step+i]
        particle_velocity[i+1], particle_position[i+1] = finding_velocity_and_position(mass_of_particle, particle_velocity[i], particle_position[i], B_field, new_time_step)

    x_array = np.zeros(len(particle_position))
    y_array = np.zeros(len(particle_position))
    z_array = np.zeros(len(particle_position))
    
    for i in range(0, len(particle_position)):
        x_array[i] = particle_position[i][0]
        y_array[i] = particle_position[i][1]
        z_array[i] = particle_position[i][2]

    return particle_velocity, particle_position, x_array, y_array, z_array

def plot_particle_motion(v_init, B_field, len_of_grid):
	particle_velocity, particle_position, x_array, y_array, z_array = equations_of_motion(v_init, B_field)

	mpl.rcParams['legend.fontsize'] = 10

	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.plot(x_array, y_array, z_array, label='particle trajectory')
        ax.set_xlim(0, length_of_grid)
        ax.set_ylim(0, length_of_grid)
        ax.set_zlim(0, length_of_grid)
	ax.set_xlabel('x position')
	ax.set_ylabel('y position')
	ax.set_zlabel('z position')
	ax.legend()
	plt.show()
