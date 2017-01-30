"""
DOCUMENTATION

PURPOSE
-------
	Calculates the characteristic length scale of the inputed grid
	with the inputted parameter.

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

def test_cases(grid, inputted_parameter, roof=1e22):
	"""
	Function takes the derivative for each value, checks
	if it equals 0 then converts values into the
	characteristic length scale
	"""
	ngrid = len(grid)
	characteristic_length_scale = []
	for i in range(0, ngrid-1):
		dq = inputted_parameter[0+i] - inputted_parameter[1+i]
		dx = grid[0+i] - grid[1+i]
		dq_dx = dq/dx
		"""
		if dq_dx != 0:
			grid_array = np.zeros(ngrid)
			characteristic_length_scale.append([x + roof for x in grid_array])
		"""
		if dq_dx != 0:
			characteristic_length_scale.append(np.abs((dq_dx/inputted_parameter[i])**-1))
		if dq_dx == 0:
			dq_dx = 1e-30
			characteristic_length_scale.append(np.abs((dq_dx/inputted_parameter[i])**-1))
		if inputted_parameter[i] == 0:
			inputted_parameter[i] = 1e-30
			characteristic_length_scale.append(np.abs((dq_dx/inputted_parameter[i])**-1))
	return np.asarray(characteristic_length_scale)

def get_loop_cls(grid, inputted_parameter, roof=1e22):
	"""
	Function checks if the length scale is finite then outputs the final characteristic
	length scales
	"""
	characteristic_length_scale = test_cases(grid, inputted_parameter, roof=1e22)
	length_of_cls = len(characteristic_length_scale)
	for i in range(0, length_of_cls):
		if np.isfinite(characteristic_length_scale[i]) == False:
			characteristic_length_scale[i] = roof
	return characteristic_length_scale 
