
import numpy as np
import statsmodels.api as sm

np.random.seed(0)

class Fit:

	def __init__(self, xx, yy, n_breakpoints, start_values):

		self.xx = xx
		self.yy = yy
		self.n_breakpoints = n_breakpoints
		self.start_values = start_values

		self.max_iterations = 10
		self.tolerance=0.1

		self.breakpoint_history = [start_values]
		self.params_history = []
		self.stop = False

		self.fit()


	def _is_valid_xx(self, xx):
		return xx

	def _is_valid_yy(self, yy):
		return yy

	def _is_valid_n_breakpoints(self, n_breakpoints):
		return n_breakpoints

	def _is_valid_start_values(self, start_values):
		return start_values

	def fit(self):
		while not self.stop:
			next_breakpoints, params = breakpoint_fit(self.xx, self.yy, self.breakpoint_history[-1]) 
			self.breakpoint_history.append(next_breakpoints)
			self.params_history.append(params)
			self.stop_or_not()

	def stop_or_not(self):

		# Stop if maximum iterations reached
		if len(self.breakpoint_history)>self.max_iterations:
			self.stop=True
		# Stop if the last change was small
		# How to do this for multiple breakpoints

"""
The breakpoint fit function is seperate to the main class
Easier for testing and re-use etc. 
Transofrms the data based on step functions
Fits based on Muggeo's method 
"""


def breakpoint_fit(xx, yy, current_breakpoints):
	"""
	Fit the linear approximation given the current breakpoint guesses
	Return the next breakpoints and the params from the fit
	The params are of the form [a, c, beta_hats, gamma_hats]
	"""

	print(current_breakpoints)

	Z = np.array([xx])
	# Convert data based on breakpoints
	UU = [(xx - bp) * np.heaviside(xx- bp, 1) for bp in current_breakpoints]
	VV = [np.heaviside(xx- bp, 1) for bp in current_breakpoints]

	print(Z, UU, VV)

	Z = np.concatenate((Z, UU, VV))
	Z = Z.T	
	Z = sm.add_constant(Z, has_constant='add')

	results = sm.OLS(endog=yy, exog=Z).fit()
	
	# First two params are a and c in the line equation
	# Beta hats are the next group of params, same length as the number of breakpoints
	beta_hats = results.params[2:2+len(current_breakpoints)]
	# Gamma hats are the last group of params, same length as the number of breakpoints
	gamma_hats = results.params[2+len(current_breakpoints):]
	# The next breakpoints are calculated iteratively
	next_breakpoints = current_breakpoints - gamma_hats/beta_hats

	print(next_breakpoints)

	return next_breakpoints, results.params








def breakpoint_fit_1_bp(xx, yy, current_breakpoints):

	current_breakpoint_1 = current_breakpoints[0]

	# Prepare data for OLS regression in SM
	A = xx
	B = (xx - current_breakpoint_1) * np.heaviside(xx- current_breakpoint_1, 1) 
	C = np.heaviside(xx - current_breakpoint_1, 1)

	print("A:" , A)

	print("B: " , B)



	Z = np.array([A, B , C]).T

	# By default add_constant won't add a column of 1s if it already exists, which it sometimes does
	# if the breakpoint is out of the data range. 
	# Change behaviour using "has_constant='add'" so we get enough columns
	Z = sm.add_constant(Z, has_constant='add')

	# Fit with SM OLS
	results = sm.OLS(endog=yy, exog=Z).fit()

	# Extract the coefficient estimates and use to get the next breakpoint
	beta_hat = results.params[2]
	print("Beta hat is ", beta_hat)
	gamma_hat = results.params[3]
	next_breakpoint_1 = current_breakpoint_1 - gamma_hat/beta_hat
	
	print(next_breakpoint_1, results.params)

	return next_breakpoint_1, results.params


def test_on_data_1():

	alpha = -4
	beta_1 = -2
	intercept = 100
	breakpoint_1 = 7

	n_points = 200

	xx = np.linspace(0, 20, n_points)

	yy = intercept + alpha*xx + beta_1 * np.maximum(xx - breakpoint_1, 0) + np.random.normal(size=n_points)

	

	bp_fit = Fit(xx, yy, n_breakpoints=1, start_values=[5])


	print(bp_fit.breakpoint_history)







def test_on_data_2():

	filename = "data/test_data_simple_1_bp.csv"

	# Results from running r segmented pacakge on the same data
	bp_result = 39
	intercept = 3.19471 
	alpha_1 = -0.08223 
	alpha_2 = 0.08410

	import pandas as pd
	df = pd.read_csv(filename, header=0)
	xx = df["x"].to_numpy()
	yy = df["y"].to_numpy()
	
	#next_bp, params = breakpoint_fit_1_bp(xx,yy,[35.122])
	

	bp_fit = Fit(xx, yy, n_breakpoints=1, start_values=[25])


	print(bp_fit.breakpoint_history)	

test_on_data_1()