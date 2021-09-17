Piecewise (aka segmented) regression in Python. Simultaneously find breakpoints and straightline segments between those breakpoints. Based on Muggeo "Estimating regression models with unknown break-points" (2003)


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

Here is a more detailed example. We start off genreating some data with a breakpoint, for demonstration purposes: ::

	import piecewise_regression
	import numpy as np

	alpha_1 = -4    
	alpha_2 = -2
	intercept = 100
	breakpoint_1 = 7
	n_points = 200
	np.random.seed(0)
	xx = np.linspace(0, 20, n_points)
	yy = intercept + alpha_1*xx + (alpha_2-alpha_1) * np.maximum(xx - breakpoint_1, 0) + np.random.normal(size=n_points)


Now we fit the model: ::

    # Given some data, fit the model
    bp_fit = Fit(xx, yy, start_values=[5], n_breakpoints=1)

    # Print a summary of the fit
    bp_fit.summary()

Example output: ::

	                    Breakpoint Regression Results                     
	====================================================================================================
	No. Observations                      200
	No. Model Parameters                    4
	Degrees of Freedom                    196
	Res. Sum of Squares               193.264
	Total Sum of Squares              46201.8
	R Squared                        0.995817
	Adjusted R Squared               0.995731
	Converged:                           True
	====================================================================================================
	====================================================================================================
	                    Estimate      Std Err            t        P>|t|       [0.025       0.975]
	----------------------------------------------------------------------------------------------------
	const                100.726        0.244       413.63     3.1e-290       100.25       101.21
	alpha1              -4.21998       0.0653      -64.605    4.37e-134      -4.3488      -4.0912
	beta1                2.18914       0.0689       31.788            -       2.0533        2.325
	breakpoint1          6.48706        0.137            -            -       6.2168       6.7573
	----------------------------------------------------------------------------------------------------
	These alphas(gradients of segments) are estimated from betas(change in gradient)
	----------------------------------------------------------------------------------------------------
	alpha2              -2.03084       0.0218      -93.068    3.66e-164      -2.0739      -1.9878
	====================================================================================================

	Davies test for existence of at least 1 breakpoint: p=5.13032e-295 (e.g. p<0.05 means reject null hypothesis of no breakpoints at 5% significance)

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

.. image:: https://raw.githubusercontent.com/chasmani/piecewise-regression/master/paper/example2.png
    :alt: fit-example-plot

How It Works
======================

The package implements Muggeo's iterative algorithm (Muggeo "Estimating regression models with unknown break-points" (2003)), to quickly find breakpoints. The Fit method also implements a non-parametric bootstrap restarting to escape local minima, this can be controlled with `n_boot`. To run the Fit without bootstrap restarting, set `n_boot=0`. Muggeo's algorthm does not always converge. In this case, the Fit method will keep trying to find a fit using bootstrap restarting `n_boot` times. 

If you don't have good guesses for inital breakpoints, you can just set the number of e.g. `n_breakpoints=3`. in this case the algorithm will randomly generate starting breakpoints until it finds a slution that converges (up to `n_boot` times). This is a good option if the algorithm is otherwise not converging. 

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

There are also a series of simluation tests that check the estimates have realistic confidence intervals, and the Davies test gives realistic p-values. These can be found in the folder "tests"

Documentation
==============
`Full docs, including an API reference. <https://piecewise-regression.readthedocs.io/en/latest/>`_
