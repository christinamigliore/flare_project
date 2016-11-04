"""
The Monte Carlo approach to approximating the value of Pi.

This method involves using indirect sampling of randomly generated (x, y) points to 
calculate the ratio between the points that fall into a square of area 4 and a circle 
with radius 1 (with the assumption the points fall within the square). This ratio should 
equal pi/4, from this pi can be approximated.

Parameters:

	tries = number of indirect samples

Returns:

	The approximation of pi and the percent error
"""

import random
import math
import numpy as np

def approximate_pi(tries=10**6):
    success = 0
    for _ in range(tries):
        x, y = random.random(), random.random()
        if x ** 2 + y ** 2 < 1:
            success += 1
    return float(4 * success) / tries
    
if __name__ == '__main__':
	approximate_pi()
	#print "% error = "+np.str(np.abs((math.pi - approximate_pi())*100/math.pi))
    
