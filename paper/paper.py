from piecewise_regression import Fit, ModelSelection
import numpy as np
import matplotlib.pyplot as plt

import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))


def plot_basic_example():
    """
    Example for some data

    """
    np.random.seed(1)

    alpha = 4
    beta_1 = -8
    beta_2 = -2
    beta_3 = 5
    intercept = 100
    breakpoint_1 = 5
    breakpoint_2 = 11
    breakpoint_3 = 16

    n_points = 200
    noise = 3

    xx = np.linspace(0, 20, n_points)

    yy = intercept + alpha*xx
    yy += beta_1 * np.maximum(xx - breakpoint_1, 0)
    yy += beta_2 * np.maximum(xx - breakpoint_2, 0)
    yy += beta_3 * np.maximum(xx - breakpoint_3, 0)
    yy += np.random.normal(size=n_points) * noise

    bp_fit = Fit(xx, yy, start_values=[3, 7, 10])

    bp_fit.summary()

    bp_fit.plot_data(color="grey", s=20)
    bp_fit.plot_fit(color="red", linewidth=4)
    bp_fit.plot_breakpoints()
    bp_fit.plot_breakpoint_confidence_intervals()

    plt.xlabel("x")
    plt.ylabel("y")

    plt.savefig("example.png", dpi=300)

    plt.show()


def plot_basic_example_2():

    # Generate some test data with 1 breakpoint
    alpha_1 = -4
    alpha_2 = -2
    intercept = 100
    breakpoint_1 = 7
    n_points = 200
    np.random.seed(0)

    xx = np.linspace(0, 20, n_points)
    yy = intercept + alpha_1*xx + \
        (alpha_2-alpha_1) * np.maximum(xx - breakpoint_1, 0) + \
        np.random.normal(size=n_points)

    # Given some data, fit the model
    bp_fit = Fit(xx, yy, start_values=[5], n_breakpoints=1)

    # Print a summary of the fit
    bp_fit.summary()

    # Plot the data, fit, breakpoints and confidence intervals
    bp_fit.plot_data(color="grey", s=20)
    # Pass in standard matplotlib keywords to control any of the plots
    bp_fit.plot_fit(color="red", linewidth=4)
    bp_fit.plot_breakpoints()
    bp_fit.plot_breakpoint_confidence_intervals()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig("example2.png")
    plt.show()
    plt.close()


def model_selection_basic_example():

    # Generate some test data with 1 breakpoint
    alpha_1 = -4
    alpha_2 = -2
    intercept = 100
    breakpoint_1 = 7
    n_points = 200
    np.random.seed(0)

    xx = np.linspace(0, 20, n_points)
    yy = intercept + alpha_1*xx + \
        (alpha_2-alpha_1) * np.maximum(xx - breakpoint_1, 0) + \
        np.random.normal(size=n_points)

    # Given some data, fit the model
    ms = ModelSelection(xx, yy, max_breakpoints=6)
    print(ms)


if __name__ == "__main__":
    plot_basic_example()
