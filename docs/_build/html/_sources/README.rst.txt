Piecewise (aka segmented) regression in Python. Simultaneously find breakpoints and fit straightline segments between those breakpoints. Based on Muggeo "Estimating regression models with unknown break-points" (2003)


Installation
========================

You can install piecewise-regression from `PyPI <https://pypi.org/project/piecewise-regression/>`_

    pip install piecewise-regression

The package was developed and tested on Python 3.7.

Getting started
========================

The package requires some x and y data to fit. You also need to specify either a) some initial breakpoint guesses as `start_values` or b) how many breakpoints you want to fit as `n_breakpoints` (or both). Here is a very simple example: ::

	import piecewise_regression
	pw_fit = piecewise_regression.Fit(x, y, n_breakpoints=2)
	pw_fit.summary()

Example
========================

Here is a more detailed example. We start off generating some data with 3 breakpoints, for demonstration purposes: ::

	import piecewise_regression
	import numpy as np

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
	noise=3

	xx = np.linspace(0, 20, n_points)

	yy = intercept + alpha*xx 
	yy += beta_1 * np.maximum(xx - breakpoint_1, 0) 
	yy += beta_2 * np.maximum(xx - breakpoint_2, 0)  
	yy += beta_3 * np.maximum(xx - breakpoint_3, 0)
	yy += np.random.normal(size=n_points) * noise


Now we fit the model: ::

	# Given some data, fit the model
	bp_fit = piecewise_regression.Fit(xx, yy, start_values=[3,7,10])

	# Print a summary of the fit
	bp_fit.summary()

Example output: ::

	                    Breakpoint Regression Results                     
	====================================================================================================
	No. Observations                      200
	No. Model Parameters                    8
	Degrees of Freedom                    192
	Res. Sum of Squares               1448.83
	Total Sum of Squares              77195.4
	R Squared                        0.981232
	Adjusted R Squared               0.980446
	Converged:                           True
	====================================================================================================
	====================================================================================================
	                    Estimate      Std Err            t        P>|t|       [0.025       0.975]
	----------------------------------------------------------------------------------------------------
	const                99.3134        0.765       129.77    5.73e-189       97.804       100.82
	alpha1               4.24777        0.268       15.861      9.3e-37       3.7196        4.776
	beta1               -8.26555        0.347      -23.848            -      -8.9492      -7.5819
	beta2               -1.80202        0.325      -5.5523            -      -2.4422      -1.1619
	beta3                5.21108        0.456       11.423            -       4.3113       6.1109
	breakpoint1          4.99612        0.129            -            -       4.7419       5.2503
	breakpoint2           10.573        0.581            -            -        9.428       11.718
	breakpoint3          16.0829        0.223            -            -       15.644       16.522
	----------------------------------------------------------------------------------------------------
	These alphas(gradients of segments) are estimated from betas(change in gradient)
	----------------------------------------------------------------------------------------------------
	alpha2              -4.01778         0.22      -18.262     7.58e-44      -4.4517      -3.5838
	alpha3              -5.81981        0.239      -24.391     1.02e-60      -6.2904      -5.3492
	alpha4             -0.608729        0.389      -1.5656        0.119      -1.3756      0.15816
	====================================================================================================

	Davies test for existence of at least 1 breakpoint: p=0.0 (e.g. p<0.05 means reject null hypothesis of no breakpoints at 5% significance)

There are also tools for plotting data: ::

	import matplotlib.pyplot as plt

	# Plot the data, fit, breakpoints and confidence intervals
	bp_fit.plot_data(color="grey", s=20)
	# Pass in standard matplotlib keywords to control any of the plots
	bp_fit.plot_fit(color="red", linewidth=4) 
	bp_fit.plot_breakpoints()
	bp_fit.plot_breakpoint_confidence_intervals()
	plt.xlabel("x")
	plt.ylabel("y")
	plt.show()
	plt.close()

.. image:: ../paper/example.png
    :alt: example-plot

How It Works
======================

The package implements Muggeo's iterative algorithm (Muggeo "Estimating regression models with unknown break-points" (2003)), to quickly find breakpoints. 

This iteartive method does not always converge to a global optimal solution, and can instead converge to a local optima or not converge at all. For this reason the Fit method also implements a non-parametric bootstrap restarting to escape local minima, this can be controlled with `n_boot`. To run the Fit without bootstrap restarting, set `n_boot=0`. If Muggeo's algorthm has not converged, the Fit method will keep trying to find a fit using bootstrap restarting `n_boot` times. 

If you don't have good guesses for inital breakpoints, you can just set the number of e.g. `n_breakpoints=3`. in this case the algorithm will randomly generate starting breakpoints until it finds a solution that converges (up to `n_boot` times). This is a good option if the algorithm is otherwise not converging. 

Model Selection
==========================

in addition to the main Fit tool, the package also offers a `ModelSelection` option based on the Bayesian Information Criterion. This is experimental and not as thorough as the main Fit tool: ::

	ms = ModelSelection(x, y, max_breakpoints=6)

This gives the following example output: ::

	                 Breakpoint Model Comparision Results                 
	====================================================================================================
	n_breakpoints            BIC    converged          RSS 
	----------------------------------------------------------------------------------------------------
	0                     421.09         True       1557.4 
	1                     14.342         True       193.26 
	2                     22.825         True       191.23 
	3                     24.169         True       182.59 
	4                     29.374         True       177.73 
	5                                   False              
	6                                   False              

	Minimum BIC (Bayesian Information Criterion) suggests the best model 



Testing
============

The package includes comprehensive tests.

To run all tests, from the main directory run: ::
	
	python3 -m "nose"

Note: This requires nosetests, can be downloaded from apt with: ::

	sudo apt install python3-nose

There are also a series of simulation tests that check the estimates have realistic confidence intervals, and the Davies test gives realistic p-values. These can be found in the folder "tests".

Documentation
==============
`Full docs, including an API reference. <https://piecewise-regression.readthedocs.io/en/latest/>`_
