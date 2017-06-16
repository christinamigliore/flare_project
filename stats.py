"""
DOCUMENTATION

PURPOSE
-------
 	Program randomly samples a uniform distribution between 0
 	and 1 then alignes these points with the y-axis to create a
 	line parallel to the x axis. The intersection between this
 	line and a lognormal cdf is found and the (x, y) coordinates
 	are saved. These x coordinates are then plotted on a
 	lognormal pdf as a histogram in order to populate the pdf.

INPUTS
------
	num_samples: integer
		Number of randomly sampled points to find.

        
OUTPUTS
-------
	Plot of lognormal pdf with the randomly found
	points as a histogram.

"""
from scipy.stats import lognorm
from scipy.special import erfc, erfcinv
import matplotlib.pyplot as plt
import math
import numpy as np
import seaborn as sns
sns.set(color_codes=True)

#Data Points
sigma = 1
x = np.linspace(lognorm.ppf(0.01, sigma),lognorm.ppf(0.99, sigma), 1000)
lognorm_cdf = np.zeros(len(x))
lognorm_pdf = np.zeros(len(x))
for i in range(0, len(x)):
    lognorm_cdf[i] = (0.5 * erfc(-math.log(x[i])/(sigma*math.sqrt(2))))
    lognorm_pdf[i] = (1 / (sigma*x[i]*math.sqrt(2*math.pi)) * math.exp(-1/2*(math.log(x[i])/sigma)**2))

def sample_cdf(num_samples):
    """
    Function randomly samples a uniform distribution between 0 and 1 based
    on the inputted num_samples. For each random values, an array is created
    that is made up of this values and the intersection index between this array 
    and the lognormal plot is found. An array for the x and y intersection points 
    are created. Returns the list of the (x, y) intersection coordinates. 
    """
    y_vals = np.random.uniform(0, 1, num_samples)
    x_vals = np.zeros(len(y_vals))
    for i in range(0, num_samples):
        x_vals[i] = math.exp(sigma*math.sqrt(2)*erfcinv(2*y_vals[i]))
    return x_vals

def plot_pdf():
    x_vals = sample_cdf(1000)
    sns.distplot(x_vals, bins = 150)
    plt.xlim(0, 15)
    plt.xlabel("x")
    plt.ylabel("f (x)")
    plt.title("Lognormal PDF with Simulated Data")
    plt.show()
