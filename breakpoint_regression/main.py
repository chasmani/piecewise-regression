
import math

import numpy as np
import statsmodels.api as sm
import scipy.stats
import matplotlib.pyplot as plt



try:
	import breakpoint_regression.davies as davies 
	import breakpoint_regression.r_squared_calc as r_squared_calc
except:
	import davies
	import r_squared_calc


class Fit:

	def __init__(self, xx, yy, start_values, max_iterations=30, tolerance=10**-5,
		min_distance_between_breakpoints=0.02, min_distance_to_edge=0.02, verbose=True):

		self.xx = self._validate_xx(xx)
		self.yy = self._validate_yy(yy)

		self.max_iterations = max_iterations	
		self.tolerance = tolerance
		# In terms of proportion of the data range
		self.min_distance_between_breakpoints = min_distance_between_breakpoints
		# In terms of quantiles
		self.min_distance_to_edge = min_distance_to_edge
		self.verbose = verbose

		self.start_values = self._validate_breakpoint_values(start_values, is_start_values=True)

		self.n_breakpoints = len(start_values)		

		self.breakpoint_history = [start_values]
		# In the format [c, alpha, betas. gammas]
		self.params_history = []
		self.covariance_history = []

		# Constant estimate
		self.const_estimate = None
		self.const_standard_error = None
		self.const_confidence_interval = None

		# Alphas are slopes of the line segments
		self.alpha_estimates = None
		self.alpha_standard_errors = None
		self.alpha_confidence_intervals = None

		# Betas are differences between slopes in line segments
		self.beta_estimates = None
		self.beta_standard_errors = None
		self.beta_confidence_intervals = None

		# Breakpoint estimates
		self.breakpoint_estimates = None
		self.breakpoint_standard_errors = None
		self.breakpoint_confidence_intervals = None

		# Davies p-value is approximately the probability of the data given there is no breakpoint
		self.davies = None

		# R squared
		self.r_squared = None
		self.adjusted_r_squared = None

		self.stop = False

		if self.verbose:
			print("Input data seems okay")

		self.fit()


	def _validate_xx(self, xx):
		"""
		Allowed types:
			List of integers of floats
			Numpy array of integers or floats
		"""
		# If its a list, convert it to a numpy array
		if isinstance(xx, list):
			xx = np.array(xx)

		if isinstance(xx, np.ndarray):
			pass

		return xx

	def _validate_yy(self, yy):
		return yy


	def _validate_breakpoint_values(self, breakpoint_values, is_start_values=False):

		breakpoint_values = self._validate_breakpoint_values_within_range(breakpoint_values, is_start_values=is_start_values)
		breakpoint_values = self._validate_breakpoint_values_far_apart(breakpoint_values, is_start_values=is_start_values)
		return breakpoint_values

	def _validate_breakpoint_values_within_range(self, breakpoints, is_start_values=False):

		min_allowed_bp = np.quantile(self.xx, self.min_distance_to_edge)
		max_allowed_bp = np.quantile(self.xx, 1-self.min_distance_to_edge)

		if is_start_values:
			value_error_text = """
				Invalid start guesses for the breakpoints
				"""
		else:
			value_error_text = """
				During the algorithm, the breakpoint values became invalid. 
				This suggests that the algorithm is not converging on good breakpoints. 
				"""

		value_error_text = value_error_text + """
			Breakpoint values are outside the allowed range of {} to {}.
			The allowed range can be changed using min_distance_to_edge, 
			min_distance_to_edge is the allowed distance to the edge in terms of quantiles of x.
			Try changing the initial breakpoint guesses, or use less breakpoints.
			The initial guesses should not be too close together, or too close to the edge of the data. 
			""".format(min_allowed_bp, max_allowed_bp)

		# Breakpoints have to be within the range of the data, plus or minus some distance
		for bp in breakpoints:
			if bp <= min_allowed_bp or bp >= max_allowed_bp:
				raise ValueError(value_error_text)

		return breakpoints

	def _validate_breakpoint_values_far_apart(self, breakpoints, is_start_values=False):

		if is_start_values:
			value_error_text = """
				Invalid start guesses for the breakpoints
				"""
		else:
			value_error_text = """
				During the algorithm, the breakpoint values became invalid. 
				This suggests that the algorithm is not converging on good breakpoints. 
				"""

		value_error_text = value_error_text + """
			Breakpoint values are too close together.
			The allowed distance can be changed using min_distance_between_breakpoints
			min_distance_between_breakpoints is the allowed distance in terms of a proportion of the range of x.
			Try changing the initial breakpoint guesses, or use less breakpoints.
			The initial guesses should not be too close together, or too close to the edge of the data. 
			"""

		# If the breakpoints are too close together, stop the algorithm and raise an error
		min_distance = np.diff(np.sort(breakpoints))
		if min_distance <= self.min_distance_between_breakpoints * np.ptp(self.xx):
			raise ValueError(value_error_text)

		return breakpoints


	def fit(self):
		# Run the breakpoint iterative procedure
		if self.verbose:
			print("Running algorithm . . . ")
		while not self.stop:
			next_breakpoints, params, cov = self.breakpoint_fit(self.xx, self.yy, self.breakpoint_history[-1]) 
			if self.verbose:
				print("Current breakpoints are {}".format(next_breakpoints))
			# Check new breakpoints are valid
			next_breakpoints = self._validate_breakpoint_values(next_breakpoints)

			# Save everything
			self.breakpoint_history.append(next_breakpoints)
			self.params_history.append(params)
			self.covariance_history.append(cov)
			self.stop_or_not()

		# Get final values, stanrd errors and confidence intervals for variables of interest
		# intercept, line segment gradients, gradient differences and breakpoints
		self.calculate_all_estimates()
		self.calculate_all_standard_errors()
		self.calculate_all_confidence_intervals()
		self.davies = davies.davies_test(self.xx, self.yy)
		self.calculate_r_squared()

	def stop_or_not(self):

		# Stop if maximum iterations reached
		if len(self.breakpoint_history)>self.max_iterations:
			if self.verbose:
				print("Max iterations reached. Stopping.")
			self.stop=True

		# Stop if tolerance reached - small change between this and last breakpoints
		breakpoint_differences = self.breakpoint_history[-2] - self.breakpoint_history[-1]
		if np.max(np.abs(breakpoint_differences)) <= self.tolerance:
			if self.verbose:
				print("Algorithm has converged on breakpoint values. Stopping.")
			self.stop=True

		# Stop if the algorithm is iterating back to previous values, within tolerance
		if len(self.breakpoint_history) > 2:
			breakpoint_two_step_differences = self.breakpoint_history[-3] - self.breakpoint_history[-1]
			if np.max(np.abs(breakpoint_two_step_differences)) <= self.tolerance:
				# Take an average of the last two values
				breakpoints = (self.breakpoint_history[-2] + self.breakpoint_history[-1])/2
				self.breakpoint_history.append(breakpoints)
				self.stop=True
				if self.verbose:
					print("Algorithm is iterating between breakpoint values. Stopping.")
					print("Final breakpoints are ", self.breakpoint_history[-1])


		# Stop if the last change was small
		# How to do this for multiple breakpoints

	
	def breakpoint_fit(self, xx, yy, current_breakpoints):
		"""
		Fit the linear approximation given the current breakpoint guesses
		Return the next breakpoints and the params from the fit
		The params are of the form [c, a, beta_hats, gamma_hats]
		Keep this as a function without referencing self.xx etc, for easier testing
		"""


		Z = np.array([xx])
		# Convert data based on breakpoints
		UU = [(xx - bp) * np.heaviside(xx- bp, 1) for bp in current_breakpoints]
		VV = [np.heaviside(xx- bp, 1) for bp in current_breakpoints]

		
		Z = np.concatenate((Z, UU, VV))
		Z = Z.T	
		Z = sm.add_constant(Z, has_constant='add')

		results = sm.OLS(endog=yy, exog=Z).fit()

		print(results.summary())

		cov = results.cov_params()
	
		# First two params are a and c in the line equation
		# Beta hats are the next group of params, same length as the number of breakpoints
		beta_hats = results.params[2:2+len(current_breakpoints)]
		# Gamma hats are the last group of params, same length as the number of breakpoints
		gamma_hats = results.params[2+len(current_breakpoints):]
		# The next breakpoints are calculated iteratively
		next_breakpoints = current_breakpoints - gamma_hats/beta_hats

		unique_sorted_xx = np.sort(np.unique(xx))

		return next_breakpoints, results.params, cov


	def calculate_all_estimates(self):
		"""
		Save all params as object variables
		"""
		params = self.params_history[-1]
		
		# Can get most of the variables directly
		self.const_estimate = params[0]
		self.beta_estimates = params[2:self.n_breakpoints+2]		
		self.breakpoint_estimates = self.breakpoint_history[-1]

		# The slopes need to be calculated by adding up the betas and first alpha
		alphas = []
		for alpha_n in range(self.n_breakpoints + 1):
			alpha = np.sum(params[1:alpha_n+2])
			alphas.append(alpha)
		self.alpha_estimates = alphas


	def get_alpha_standard_errors(self):

		cov_matrix = self.covariance_history[-1]

		# Alphas are calculated as the sum of alpha_0 and betas up that part of the regression line 
		# The var of each alpha is the sum of the covariance matrix up to it
		# Removing the intercept column and row. 
		# var(alpha_k) = var(alpha_1) + sum var(beta_j) + 2*sum_{i=1,j=2}^k *cov(alpha, betas))
		# var(alpha_k) = sum_{i, j} cov(alpha and betas)
		alpha_vars = []
		for alpha_n in range(self.n_breakpoints + 1):
			alpha_cov_matrix = cov_matrix[1:alpha_n+2, 1:alpha_n+2]	
			alpha_vars.append(np.sum(alpha_cov_matrix))

		alpha_ses = np.sqrt(alpha_vars)	
		return alpha_ses	


	def get_bp_standard_errors(self):

		# bp = gamma/beta + bp_0
		# Variance of bp estimator found using ratio/delta method
		# See the accompanying paper for clarification

		cov_matrix = self.covariance_history[-1]
		params = self.params_history[-1]

		bp_vars = []

		# For each breakpoint, calcaulte the variance of the estimator
		for bp_n in range(self.n_breakpoints):
			beta_index = 2 + bp_n
			gamma_index = 2 + self.n_breakpoints + bp_n

			beta = params[beta_index]
			gamma = params[gamma_index]
			gamma_var = cov_matrix[gamma_index, gamma_index]
			beta_var = cov_matrix[beta_index, beta_index]
			gamma_beta_covar = cov_matrix[beta_index, gamma_index]

			# Differs from Muggeo (2003) in the minus sign before the covar term - I think Muggeo had a mistake
			bp_var = (gamma_var + beta_var * (gamma/beta)**2 - 2*(gamma/beta)*gamma_beta_covar)/(beta**2)
			bp_vars.append(bp_var)

		bp_ses = np.sqrt(bp_vars)
		return bp_ses


	def calculate_all_standard_errors(self):
		"""
		Calcaulte standrd errors for all the variables of interest
		"""
		# Covariance matrix is [c, alpha, betas, gammas]
		# Constant variance is jsut the top left cell in covariance matrix
		cov_matrix = self.covariance_history[-1]
		c_var = cov_matrix[0][0]
		self.const_standard_error = np.sqrt(c_var)
		
		# Beta variances are along the diagonal of the covariance matrix
		cov_diagonal = np.diagonal(cov_matrix)
		beta_vars = cov_diagonal[2:self.n_breakpoints+2]
		self.beta_standard_errors = np.sqrt(beta_vars)

		self.alpha_standard_errors = self.get_alpha_standard_errors()
		self.breakpoint_standard_errors = self.get_bp_standard_errors()


	def calculate_all_confidence_intervals(self):
		"""
		Final confidence intervals are for [x, alphas, bps]
		"""	
		# Calcualte the confidence intervals based on Student's t distribution
		dof = len(self.xx) - 2 - 2*self.n_breakpoints
		t_const = scipy.stats.t.ppf(0.975, dof)

		self.const_confidence_interval = (self.const_estimate - t_const*self.const_standard_error, 
			self.const_estimate + t_const*self.const_standard_error)

		# For the params with maybe more than one value, use a zip to handle the list unpacking
		self.alpha_confidence_intervals = tuple(zip(self.alpha_estimates - t_const*self.alpha_standard_errors, 
			self.alpha_estimates + t_const*self.alpha_standard_errors))

		self.beta_confidence_intervals = tuple(zip(self.beta_estimates - t_const*self.beta_standard_errors, 
			self.beta_estimates + t_const*self.beta_standard_errors))

		self.breakpoint_confidence_intervals = tuple(zip(self.breakpoint_estimates - t_const*self.breakpoint_standard_errors, 
			self.breakpoint_estimates + t_const*self.breakpoint_standard_errors))
	
	def get_predicted_yy(self):

		final_params = self.params_history[-1]
		breakpoints = self.breakpoint_history[-1]
		# Extract what we need from params etc
		intercept_hat = final_params[0]
		alpha_hat = final_params[1]
		beta_hats = final_params[2:2+len(breakpoints)]

		yy_predicted = []

		yy_predicted = intercept_hat + alpha_hat*self.xx
		for bp_count in range(len(breakpoints)):
			yy_predicted += beta_hats[bp_count] * np.maximum(self.xx - breakpoints[bp_count], 0)
		
		print(list(self.yy), list(yy_predicted))

		return yy_predicted

	def calculate_r_squared(self):

		yy_predicted = self.get_predicted_yy()
		n_params = 2 * self.n_breakpoints + 2

		self.r_squared, self.adjusted_r_squared = r_squared_calc.get_r_squared(self.yy, yy_predicted, n_params)
		print(self.r_squared, self.adjusted_r_squared)


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
		# Params are in terms of [intercept, alpha, betas, gammas]
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

	def plot_breakpoint_confidence_intervals(self, **kwargs):
		"""
		Plot the breakpoint cis as shaded regions
		"""
		for bp_ci in self.breakpoint_confidence_intervals:
			plt.axvspan(bp_ci[0], bp_ci[1], alpha=0.1)



	def summary(self):
		print("Coming soon")
		header = "{:^50}\n".format("Breakpoint Regression Results")

		line_length=60
		double_line = "=" * line_length + "\n"
		single_line = "-" * line_length + "\n"

		n_obs = len(self.xx)
		dof = 2 + 2*self.n_breakpoints

		no_obs_text = "{:<20} {:>20}\n".format("No. Observations", n_obs)
		dof_text = "{:<20} {:>20}\n".format("Degress of Freedom", dof)
		r_2_text = "{:<20} {:>20.6f}\n".format("R Squared", self.r_squared)
		adj_r_2_text = "{:<20} {:>20.6f}\n".format("Adjusted R Squared", self.adjusted_r_squared)

		overview = double_line + no_obs_text + dof_text + r_2_text + adj_r_2_text + double_line


		print(header + overview)



"""
The breakpoint fit function is seperate to the main class
Easier for testing and re-use etc. 
Transforms the data based on step functions before fitting
Fits based on Muggeo's method 
"""



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

	print("COV: ", results.cov_params())

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


	print("p-value is ", bp_fit.davies)


	#print(bp_fit.breakpoint_history)

	#bp_fit.plot_data()
	#plt.show()


def test_on_data_1b():

	alpha = -4
	beta_1 = -4
	beta_2 = 4
	intercept = 100
	breakpoint_1 = 7
	breakpoint_2 = 12

	n_points = 200

	xx = np.linspace(0, 20, n_points)

	yy = intercept + alpha*xx + beta_1 * np.maximum(xx - breakpoint_1, 0) + beta_2 * np.maximum(xx-breakpoint_2, 0)  + np.random.normal(size=n_points)


	bp_fit = Fit(xx, yy, start_values=[5, 10])

	bp_fit.summary()

	bp_fit.plot_breakpoint_history()
	#plt.show()



	bp_fit.plot_data()
	bp_fit.plot_fit(color="red", linewidth=4)
	bp_fit.plot_breakpoints()
	bp_fit.plot_breakpoint_confidence_intervals()
	#plt.show()



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


	bp_fit = Fit(xx, yy, n_breakpoints=4, start_values=[3,5, 10,11])

	bp_fit.plot_data()
	bp_fit.plot_fit(color="red", linewidth=4)
	bp_fit.plot_breakpoints()
	bp_fit.plot_breakpoint_confidence_intervals()



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

if __name__=="__main__":

	np.random.seed(0)
	test_on_data_1b()