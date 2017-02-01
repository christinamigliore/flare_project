"""
DOCUMENTATION

PURPOSE
-------
	Calculates the characteristic length scale of the inputed grid
	with the inputted parameter, (q).

INPUTS
------
	grid: 1D array
		grid values

	inputted_parameter: 1D array
		quantity values

	roof: float, optionalg
		maximum value
		

OUTPUTS
-------
	c_ls: array
		characteristic length scale of the grid

"""
import numpy as np

def get_loop_cls(grid, inputted_parameter, roof=1e22):
    """
    Function takes the derivative for each value, checks
    if it equals 0 then converts values into the
    characteristic length scale
    """
    ngrid = len(grid)
    dq = np.diff(inputted_parameter)
    dx = np.diff(grid)
    dq_dx = dq/dx

    if len(dq_dx[np.where(dq_dx == 0)])+1 == ngrid:
        grid_array = np.zeros(ngrid)
        characteristic_length_scale = [x + roof for x in grid_array] 

    else:
        index_zero = np.argwhere(dq_dx == 0)
        index_zero_parameter = np.argwhere(inputted_parameter == 0)

        zero_length = len(index_zero)
        zero_parameter_length = len(index_zero_parameter)
    
        for i in range(0, zero_length):
            dq_dx[index_zero[i]] = 1e-30

        for j in range(0, zero_parameter_length):
            inputted_parameter[index_zero_parameter[j]] = 1e-30

        characteristic_length_scale = np.abs((dq_dx/inputted_parameter[: -1])**-1)

        index_infinite = np.argwhere(characteristic_length_scale)
        length_of_cls = len(characteristic_length_scale)

        for h in range(0, length_of_cls):
            if np.isfinite(characteristic_length_scale[h]) == False:
                characteristic_length_scale[h] = roof

    return characteristic_length_scale 
