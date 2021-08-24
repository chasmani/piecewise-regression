import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import math

import os, sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from piecewise_regression import Fit

def plot_basic_example():
    """
    Example for some data

    """
    np.random.seed(1)

    alpha = 4
    beta_1 = -8
    beta_2 = -2
    beta_3 = 3
    intercept = 100
    breakpoint_1 = 5
    breakpoint_2 = 11
    breakpoint_3 = 16

    n_points = 200
    noise=3

    xx = np.linspace(0, 20, n_points)

    yy = intercept + alpha*xx 
    yy += beta_1 * np.maximum(xx - breakpoint_1, 0) 
    yy += beta_2 * np.maximum(xx - breakpoint_2, 0)  
    yy += beta_3 * np.maximum(xx - breakpoint_3, 0)
    yy += np.random.normal(size=n_points) * noise

    print(xx, yy)


    bp_fit = Fit(xx, yy, start_values=[7])

    bp_fit.summary()

    bp_fit.plot_data(color="grey", s=20)
    bp_fit.plot_fit(color="red", linewidth=4)
    bp_fit.plot_breakpoints()
    bp_fit.plot_breakpoint_confidence_intervals()
    
    plt.xlabel("x")
    plt.ylabel("y")

    print("BIC is :", bp_fit.calculate_bayesian_information_criterion())

    plt.savefig("example.png", dpi=300)
    
    plt.show()

plot_basic_example()