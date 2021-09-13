
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import math

import os, sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from piecewise_regression.main import Muggeo, Fit


def check_p_values_fit():

    p_counts ={
        "alpha1":0,
        "beta1":0,
        "alpha2":0,
        "const":0
    }

 
    sample_size = 1000

    actuals = {
        "alpha1":1,
        "beta1":3,
        "alpha2":4,
        "const":5
    }

    for seed in range(sample_size):
        print("Working on {} of {} . . . ".format(seed, sample_size))
        np.random.seed(seed)

        xx_bp, yy_bp = generate_data(intercept=actuals["const"], alpha=actuals["alpha1"], beta_1=actuals["beta1"])

        bp_fit = Fit(xx_bp, yy_bp, [5], verbose=False)

        for estimate in ["alpha1", "beta1", "alpha2", "const"]:
            est_ci = bp_fit.best_muggeo.best_fit.estimates[estimate]["confidence_interval"]
            if est_ci[0] <= actuals[estimate] and est_ci[1]>= actuals[estimate]:
                p_counts[estimate] += 1

    print("Percentages in confidence interval:")
    for estimate in ["alpha1", "beta1", "alpha2", "const"]:        
        print("{} : {}".format(estimate, p_counts[estimate]/sample_size))
    print("This should be approximately 95%")    



def generate_data(intercept, alpha, beta_1):

    intercept = 5
    alpha = 1
    beta_1 = 3
    n_points = 50
    noise = 1
    breakpoint_1 = 2

    xx_bp = np.linspace(-9.5, 9.5, n_points)
    
    yy_bp = intercept + alpha*xx_bp + beta_1 * np.maximum(xx_bp - breakpoint_1, 0) + noise*np.random.normal(size=len(xx_bp))


    return xx_bp, yy_bp





if __name__=="__main__":
    check_p_values_fit()