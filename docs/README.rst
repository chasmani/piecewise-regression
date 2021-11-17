==========================================================
piecewise-regression (aka segmented regression) in python
==========================================================
:piecewise-regression: fitting straight line models with breakpoints
:Author: Charlie Pilgrim
:Version: 1.0.4
:Github: https://github.com/chasmani/piecewise-regression
:Documentation: https://piecewise-regression.readthedocs.io/en/master/index.html

.. image:: https://github.com/chasmani/piecewise-regression/actions/workflows/python-package.yml/badge.svg
   :target: https://github.com/chasmani/piecewise-regression/actions/workflows/python-package.yml
   :alt: Build Status
.. image:: https://codecov.io/gh/chasmani/piecewise-regression/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/chasmani/piecewise-regression
   :alt: Test Coverage Status
.. image:: https://readthedocs.org/projects/piecewise-regression/badge/?version=latest
   :target: https://piecewise-regression.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://badge.fury.io/py/piecewise-regression.svg
   :target: https://badge.fury.io/py/piecewise-regression
   :alt: PyPi location
.. image:: https://joss.theoj.org/papers/b64e5e7d746efc5d91462a51b3fc5bf8/status.svg
   :target: https://joss.theoj.org/papers/b64e5e7d746efc5d91462a51b3fc5bf8
   :alt: Review Status
.. image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://github.com/chasmani/piecewise-regresssion/blob/master/LICENSE
   :alt: License information

|

Easy-to-use piecewise regression (aka segmented regression) in Python. For fitting straight lines to data where there is one or more changes in gradient (known as breakpoints). Based on Muggeo "Estimating regression models with unknown break-points" (2003). 

For example:

.. image:: https://raw.githubusercontent.com/chasmani/piecewise-regression/master/paper/example.png
    :alt: basic-example-plot-github

There are some code examples below, and more in this |colab_link|.

.. |colab_link| raw:: html

   <a href="https://colab.research.google.com/drive/1Pwv6LqwZU8Zbl0VZH6cwOTwoRzm3CPPC#offline=true&sandboxMode=true" target="_blank">Google Colab Jupyter Notebook</a>


Installation
========================

You can install piecewise-regression using python's |pip_link|

.. |pip_link| raw:: html

   <a href="https://pypi.org/project/piecewise-regression" target="_blank">pip package index</a>

::

    pip install piecewise-regression

The package is tested on Python 3.7, 3.8 and 3.9.

Getting started
========================

The package requires some x and y data to fit. You also need to specify either a) some initial breakpoint guesses as `start_values` or b) how many breakpoints you want to fit as `n_breakpoints` (or both). Here is a very simple example, assuming we already have some data `x` and `y`: ::

	import piecewise_regression
	pw_fit = piecewise_regression.Fit(x, y, n_breakpoints=2)
	pw_fit.summary()

Example
========================

Here is a more detailed example. We start off generating some data with a breakpoint. This is for demonstration purposes, normally you will have your own data to fit: ::

	import piecewise_regression
	import numpy as np

	alpha_1 = -4    
	alpha_2 = -2
	constant = 100
	breakpoint_1 = 7
	n_points = 200
	np.random.seed(0)
	xx = np.linspace(0, 20, n_points)
	yy = constant + alpha_1*xx + (alpha_2-alpha_1) * np.maximum(xx - breakpoint_1, 0) + np.random.normal(size=n_points)


Now we fit the model: ::

    # Given some data, fit the model
    pw_fit = piecewise_regression.Fit(xx, yy, start_values=[5], n_breakpoints=1)

    # Print a summary of the fit
    pw_fit.summary()

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

This includes estimates for all the model variables, along with confidence intervals. The Davies test is a hypothesis test for the existence of at least one breakpoint, against the null hypothesis of no breakpoints.  

There are also tools for plotting data: ::

	import matplotlib.pyplot as plt

	# Plot the data, fit, breakpoints and confidence intervals
	pw_fit.plot_data(color="grey", s=20)
	# Pass in standard matplotlib keywords to control any of the plots
	pw_fit.plot_fit(color="red", linewidth=4) 
	pw_fit.plot_breakpoints()
	pw_fit.plot_breakpoint_confidence_intervals()
	plt.xlabel("x")
	plt.ylabel("y")
	plt.show()
	plt.close()

.. image:: https://raw.githubusercontent.com/chasmani/piecewise-regression/master/paper/example2.png
    :alt: fit-example-plot-github


You can extract data as well: ::

	# Get the key results of the fit 
	pw_results = pw_fit.get_results()
	pw_estimates = pw_results["estimates"]


How It Works
======================

The package implements Muggeo's iterative algorithm (Muggeo "Estimating regression models with unknown break-points" (2003)), to quickly find breakpoints. That method simultaneously fits breakpoint positions and the linear models for the different segments of the fit. This method is quick and it gives confidence intervals for all the model estimates. See the accompanying paper for more details.

Muggeo's method doesn't always converge on the best solution - sometimes it finds a locally optimal solution or doesn't converge at all. For this reason the Fit method also implements a process called bootstrap restarting. This involves taking a bootstrap resample of the data, then using this bootstrapped data to try and find a better solution. The number of times this runs can be controlled with `n_boot`. To run the Fit without bootstrap restarting, set `n_boot=0`.  

If you don't have good guesses for inital breakpoints, you can just set the number of e.g. `n_breakpoints=3`. in this case the algorithm will randomly generate start_values for breakpoints until it finds a solution that converges (up to `n_boot` times). This is a good option if the algorithm is otherwise not converging. Be aware that the start_values can influence the final converged model, so setting them randomly in this way may give different results on different runs, epecially if the breakpoint positions are not very clear from the data. 

Model Selection
==========================

In addition to the main Fit tool, the package also offers a `ModelSelection` option based on the Bayesian Information Criterion. This is experimental and is not as thorough as the main Fit function. In particular, the models are generated with random start_values which can influence the model fit and give different values for the BIC. The tool can be useful for exploring posisble models, but should not at this point be used to choose the best model. ::

	ms = piecewise_regression.ModelSelection(x, y, max_breakpoints=6)

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

The data of the model fits can be accessed in ::

    ms.models 

For a robust comparision, one could run the ModelSelection tools many times and take the lowest BIC for each model. 


Testing
============

The package includes comprehensive tests.

To run all tests, from the main directory run (requires the pytest library): ::
	
	pytest

To get code coverage, run (requires pytest and pytest-cov libraries): ::

	pytest --cov=./

There are also a series of simulation tests that check the estimates have realistic confidence intervals, and the Davies test gives realistic p-values. These can be found in the folder "tests-manual". 

Requirements
=============

See requirements.txt for specific version numbers. Required packages, and their uses are:

- matplotlib for plotting.
- numpy for simple data handling and data transformations.  
- scipy for statistical tests including using t-distributions and Gaussians. 
- statsmodels for performing ordinary least squares.

Community Guidelines and Contributing
===================================================

I welcome community participation:

- Open an issue on github if you want to suggest a new feature or report a bug
- If you want to make changes yourself, that is welcome via a pull request. 
- Ideally, open an issue first before making a pull request with major changes. 

Installing From Source
===========================

To install from source: ::

	git clone https://github.com/chasmani/piecewise-regression
	cd piecewise_regression
	python3 setup.py install --user


Documentation
==============
`Full docs, including an API reference. <https://piecewise-regression.readthedocs.io/en/latest/>`_




