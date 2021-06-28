
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt

np.random.seed(0)

class Fit:

	def __init__(self, xx, yy, n_breakpoints, start_values):

		self.xx = xx
		self.yy = yy
		self.n_breakpoints = n_breakpoints
		self.start_values = start_values

		self.max_iterations = 10
		self.tolerance=0.1 # Maybe set this based on gaps between data

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

	def plot_data(self, **kwargs):
		"""
		Plot the data as a scatter plot
		Passes any kwargs to the matplotlib scatter function, e.g. color="red"
		"""
		plt.scatter(self.xx, self.yy, **kwargs)

	def plot_fit(self, **kwargs):
		"""
		Plot the fitted model
		Passes any kwargs to the matplotlib plot function, e.g. color="red"
		"""
		# Get the final results from the fitted model variables
		final_params = self.params_history[-1]
		breakpoints = self.breakpoint_history[-1]
		
		# Extract what we need from params etc
		intercept_hat = final_params[0]
		alpha_hat = final_params[1]
		beta_hats = final_params[2:2+len(breakpoints)]
		
		xx_plot = np.linspace(min(self.xx), max(self.xx), 100)

		# Build the fit plot segment by segment. Betas are defined as difference in gradient from previous section
		yy_plot = intercept_hat + alpha_hat*xx_plot
		for bp_count in range(len(breakpoints)):
			yy_plot += beta_hats[bp_count] * np.maximum(xx_plot - breakpoints[bp_count], 0)

		plt.plot(xx_plot, yy_plot, **kwargs)

	def plot_breakpoints(self, **kwargs):
		"""
		Plot the breakpoint locations
		Passes kwargs to the matplotlib function, e.g. color="red"
		"""
		for bp in self.breakpoint_history[-1]:
			plt.axvline(bp, **kwargs)

	def plot_breakpoint_history(self, **kwargs):
		"""
		Plot the history of the breakpoints as they iterate
		"""
		# Change the format of the breakpoint histories to seperate lists for each breakpoint
		breakpoint_histories = zip(*self.breakpoint_history)
		# Plot a history for each breakpoint 
		count = 0
		for bh in breakpoint_histories:
			count += 1
			plt.plot(range(1, len(bh)+1), bh, label="Breakpoint {}".format(count), **kwargs)
			plt.xlabel("Iteration")
			plt.ylabel("Breakpoint")




"""
The breakpoint fit function is seperate to the main class
Easier for testing and re-use etc. 
Transforms the data based on step functions before fitting
Fits based on Muggeo's method 
"""
def breakpoint_fit(xx, yy, current_breakpoints):
	"""
	Fit the linear approximation given the current breakpoint guesses
	Return the next breakpoints and the params from the fit
	The params are of the form [c, a, beta_hats, gamma_hats]
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

	bp_fit.plot_data()
	plt.show()


def test_on_data_1b():

	alpha = -4
	beta_1 = -2
	beta_2 = 4
	intercept = 100
	breakpoint_1 = 7
	breakpoint_2 = 12

	n_points = 200

	xx = np.linspace(0, 20, n_points)

	yy = intercept + alpha*xx + beta_1 * np.maximum(xx - breakpoint_1, 0) + beta_2 * np.maximum(xx-breakpoint_2, 0)  + np.random.normal(size=n_points)


	bp_fit = Fit(xx, yy, n_breakpoints=2, start_values=[5, 10])


	print(bp_fit.breakpoint_history)

	bp_fit.plot_data()
	bp_fit.plot_fit(color="red", linewidth=4)
	bp_fit.plot_breakpoints()
	plt.show()
	plt.close()

	bp_fit.plot_breakpoint_history()
	plt.legend()
	plt.show()



def test_on_data_1c():

	alpha = -4
	beta_1 = -2
	beta_2 = 4
	beta_3 = 1
	intercept = 100
	breakpoint_1 = 7
	breakpoint_2 = 12
	breakpoint_3 = 14

	n_points = 200

	xx = np.linspace(0, 20, n_points)

	yy = intercept + alpha*xx 
	yy += beta_1 * np.maximum(xx - breakpoint_1, 0) 
	yy += beta_2 * np.maximum(xx - breakpoint_2, 0)  
	yy += beta_3 * np.maximum(xx - breakpoint_3, 0)
	yy += np.random.normal(size=n_points)


	bp_fit = Fit(xx, yy, n_breakpoints=3, start_values=[5, 10,11])


	print(bp_fit.breakpoint_history)

	bp_fit.plot_data()
	bp_fit.plot_fit(color="red", linewidth=4)
	bp_fit.plot_breakpoints()
	plt.show()
	plt.close()

	bp_fit.plot_breakpoint_history()
	plt.legend()
	plt.show()



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

test_on_data_1c()