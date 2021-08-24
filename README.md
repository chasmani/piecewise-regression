# breakpoint-regression

Breakpoint (aka segmented) regression in Python. Simultaneously find breakpoints and straightline segments between those breakpoints. Based on Muggeo "Estimating regression models with unknown break-points" (2003)


## Installation

You can install breakpoint-regression from [PyPI](https://pypi.org/project/breakpoint-regression/):

    pip install breakpoint-regression

The package is supported on Python 3.7 and above.

## How To Use

The package requires some x and y data to fit. You also need to tell it how many breakpoints you want, and starting guesses for those breakpoints. 

The package includes tools for summarising the fitted model and plotting.

	import piecewise_regression
	import matplotlib.pyplot as plt
	import numpy as np

	# Generate some test data with 1 breakpoint
	alpha_1 = -4
	alpha_2 = 1
	intercept = 100
	breakpoint_1 = 7
	n_points = 200

	xx = np.linspace(0, 20, n_points)
	yy = intercept + alpha_1*xx + (alpha_2-alpha_1) * np.maximum(xx - breakpoint_1, 0) + np.random.normal(size=n_points)

	# Given some data, fit the model
	bp_fit = piecewise_regression.Fit(xx, yy, n_breakpoints=1, start_values=[5])

	# Print a summary of the fit
	bp_fit.summary()

	# Plot the data, fit, breakpoints and confidence intervals
	bp_fit.plot_data()
	# Pass in standard matplotlib keywords to control any of the plots
	bp_fit.plot_fit(color="red", linewidth=4) 
	bp_fit.plot_breakpoints()
	bp_fit.plot_breakpoint_confidence_intervals()
	plt.show()
	plt.close()

	# Plot the history of the breakpoints during the algorithm 
	bp_fit.plot_breakpoint_history()
	plt.show()

## Testing

From the main dierctory run 
	
	python3 -m "nose"

Note: This requires nosetests, can be downloaded from apt with

	sudo apt install python3-nose
