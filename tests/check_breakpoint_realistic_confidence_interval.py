
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import math

import os, sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from piecewise_regression.main import Muggeo, Fit

def check_p_values():

    p_count = 0
    
    sample_size = 1000

    actual_bp = 2


    for seed in range(sample_size):
        print("Working on {} of {} . . . ".format(seed, sample_size))
        np.random.seed(seed)

        xx_bp, yy_bp = generate_data(actual_bp)

        bp_fit = Muggeo(xx_bp, yy_bp, [5], verbose=False)

        bp_ci = bp_fit.best_fit.estimates["breakpoint1"]["confidence_interval"]
        if bp_ci[0] <= actual_bp and bp_ci[1]>= actual_bp:
            p_count += 1  

    print(p_count)

    print("{} of {} estimates were within the confidence interval.".format(p_count, sample_size))
    print("This should be approximately 95%")

def check_p_values_fit():

    p_count = 0
    
    sample_size = 100

    actual_bp = 2


    for seed in range(sample_size):
        print("Working on {} of {} . . . ".format(seed, sample_size))
        np.random.seed(seed)

        xx_bp, yy_bp = generate_data(actual_bp)

        bp_fit = Fit(xx_bp, yy_bp, [5], verbose=False)

        bp_ci = bp_fit.best_muggeo.best_fit.estimates["breakpoint1"]["confidence_interval"]
        if bp_ci[0] <= actual_bp and bp_ci[1]>= actual_bp:
            p_count += 1  

    print(p_count)

    print("{} of {} estimates were within the confidence interval.".format(p_count, sample_size))
    print("This should be approximately 95%")    



def generate_data(breakpoint_1):

    intercept = 5
    alpha = 1
    beta_1 = 3
    n_points = 50
    noise = 1

    xx_bp = np.linspace(-9.5, 9.5, n_points)
    
    yy_bp = intercept + alpha*xx_bp + beta_1 * np.maximum(xx_bp - breakpoint_1, 0) + noise*np.random.normal(size=len(xx_bp))


    return xx_bp, yy_bp





if __name__=="__main__":
    check_p_values_fit()