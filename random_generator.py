import numpy as np

def get_random_points(number_of_points, min_lim, max_lim, dimensions):
    """
    PURPOSE
    -------
        Produces lists of randomly generated vectors of specified
        dimension and min/max values.

    INPUTS
    ------
        num_of_particles: integer
            Number of particles to be generated

        min_lim: integer
            Minimum limit to the randomly generated value

        max_lim: integer
            Maximum limit to the randomly generated value

        dimensions: integer
            Dimension of the randomly generated vector
            
    OUTPUTS
    -------
        vector: list
            List of randomly generated vectors.
    """
    vector = []
    for i in range(0, number_of_points):
        rand_nums = np.random.uniform(min_lim, max_lim, dimensions)
        vector.append(rand_nums)
    return vector 
