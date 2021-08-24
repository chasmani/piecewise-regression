
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


class NextBreakpoints:
	"""
	One iteration of Muggeo's segmented regression algorithm
	Gets the next breakpoints. 
	Also gets interesting statistics etc
	"""

	def __init__(self, xx, yy, current_breakpoints):

		self.xx = xx
		self.yy = yy
		self.current_breakpoints = current_breakpoints
		self.n_breakpoints = len(current_breakpoints)

		self.next_breakpoints = None
		self.raw_params = None
		self.covariance_matrix = None

		self.breakpoint_fit()

		# All estimate data saved in dictionary
		self.estimates = {}

		# Don't really need to do this at this point. But it is very quick and nice to have the record of these as we go 
		self.calculate_all_estimates()
		self.calculate_all_standard_errors()
		self.calculate_all_confidence_intervals()
		self.calculate_all_t_stats()

		# R squared etc
		self.residual_sum_squares = None 
		self.total_sum_squares = None
		self.r_squared = None
		self.adjusted_r_squared = None
		self.bic = None

		self.calculate_r_squared()
		self.calculate_bayesian_information_criterion()


	def breakpoint_fit(self):
		"""
		Fit the linear approximation given the current breakpoint guesses
		Return the next breakpoints and the params from the fit
		The params are of the form [c, a, beta_hats, gamma_hats]
		Keep this as a function without referencing self.xx etc, for easier testing
		"""

		Z = np.array([self.xx])
		# Convert data based on breakpoints
		UU = [(self.xx - bp) * np.heaviside(self.xx - bp, 1) for bp in self.current_breakpoints]
		VV = [np.heaviside(self.xx- bp, 1) for bp in self.current_breakpoints]

		
		Z = np.concatenate((Z, UU, VV))
		Z = Z.T	
		Z = sm.add_constant(Z, has_constant='add')

		results = sm.OLS(endog=self.yy, exog=Z).fit()

		self.raw_params = results.params
		self.covariance_matrix = results.cov_params()

		# First two params are a and c in the line equation
		# Beta hats are the next group of params, same length as the number of breakpoints
		beta_hats = results.params[2:2+len(self.current_breakpoints)]
		# Gamma hats are the last group of params, same length as the number of breakpoints
		gamma_hats = results.params[2+len(self.current_breakpoints):]
		# The next breakpoints are calculated iteratively
		self.next_breakpoints = self.current_breakpoints - gamma_hats/beta_hats



	def calculate_all_estimates(self):
		"""
		Save all params as in estimates hash table
		"""
		# Save in estimates hash table
		params = self.raw_params
	
		const_estimate = params[0]
		beta_estimates = params[2:self.n_breakpoints+2]
		breakpoint_estimates = self.next_breakpoints
		
		self.estimates["const"] = {"estimate":const_estimate}
		for bp_i in range(self.n_breakpoints):
			self.estimates["beta{}".format(bp_i+1)] = {"estimate":beta_estimates[bp_i]}
			self.estimates["breakpoint{}".format(bp_i+1)] = {"estimate":breakpoint_estimates[bp_i]}
		# Also calculate alphas
		for alpha_i in range(self.n_breakpoints + 1):
			alpha_estimate = np.sum(params[1:alpha_i+2])
			self.estimates["alpha{}".format(alpha_i+1)] = {"estimate":alpha_estimate}


	def get_alpha_standard_errors(self):

		cov_matrix = self.covariance_matrix

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

		cov_matrix = self.covariance_matrix
		params = self.raw_params

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


	def get_const_standard_error(self):
		# Covariance matrix is [c, alpha, betas, gammas]
		# Constant variance is jsut the top left cell in covariance matrix
		cov_matrix = self.covariance_matrix
		c_var = cov_matrix[0][0]
		return np.sqrt(c_var)

	def get_beta_standard_errors(self):
		# Covariance matrix is [c, alpha, betas, gammas]
		# Beta variances are along the diagonal of the covariance matrix
		cov_matrix = self.covariance_matrix
		cov_diagonal = np.diagonal(cov_matrix)
		beta_vars = cov_diagonal[2:self.n_breakpoints+2]
		return np.sqrt(beta_vars)


	def calculate_all_standard_errors(self):
		"""
		Calculate standrd errors for all the variables of interest
		Save to the estimates dictionary
		"""
		const_ses = self.get_const_standard_error()
		self.estimates["const"]["se"] = const_ses

		beta_ses = self.get_beta_standard_errors()
		bp_ses = self.get_bp_standard_errors()
		for bp_i in range(self.n_breakpoints):
			self.estimates["beta{}".format(bp_i+1)]["se"] = beta_ses[bp_i]
			self.estimates["breakpoint{}".format(bp_i+1)]["se"] = bp_ses[bp_i]
		alpha_ses = self.get_alpha_standard_errors()
		for alpha_i in range(self.n_breakpoints + 1):
			self.estimates["alpha{}".format(alpha_i+1)]["se"] = alpha_ses[alpha_i]


	def calculate_all_confidence_intervals(self):
		"""
		Confidence intervals based on t-distribution
		"""	
		# Estimates
		dof = len(self.xx) - 2 - 2*self.n_breakpoints
		t_const = scipy.stats.t.ppf(0.975, dof)

		# Iterate over the estimate dictionary, add confidence intervals to all estaimtors 
		for estimator_name, details in self.estimates.items():
			confidence_interval = (details["estimate"] - t_const*details["se"], details["estimate"] + t_const*details["se"])
			details["confidence_interval"] = confidence_interval

	def calculate_all_t_stats(self):
		"""
		Get t stats for all the estimators
		"""

		dof = len(self.xx) - 2 - 2*self.n_breakpoints
		for estimator_name, details in self.estimates.items():
			# Breakpoint t stats don't make sense
			# Don't exist in the null model - nuisance parameter
			# H_0 isn't bp=0, it's that bp doesn't exist
			if "breakpoint" in estimator_name:
				details["t_stat"] = "-"
				details["p_t"] = "-"
			else:
				t_stat = details["estimate"]/details["se"]
				p_t = scipy.stats.t.sf(np.abs(t_stat), dof)*2
				details["t_stat"] = t_stat
				details["p_t"] = p_t		


	def get_predicted_yy(self):

		params = self.raw_params
		breakpoints = self.next_breakpoints
		# Extract what we need from params etc
		intercept_hat = params[0]
		alpha_hat = params[1]
		beta_hats = params[2:2+len(breakpoints)]

		yy_predicted = []

		yy_predicted = intercept_hat + alpha_hat*self.xx
		for bp_count in range(len(breakpoints)):
			yy_predicted += beta_hats[bp_count] * np.maximum(self.xx - breakpoints[bp_count], 0)
		
		return yy_predicted

	def calculate_r_squared(self):

		yy_predicted = self.get_predicted_yy()
		n_params = 2 * self.n_breakpoints + 2

		self.residual_sum_squares, self.total_sum_squares, self.r_squared, self.adjusted_r_squared = r_squared_calc.get_r_squared(self.yy, yy_predicted, n_params)

	def calculate_bayesian_information_criterion(self):
		"""
		Assuming normal noise
		"""
		n = len(self.xx) # No. data points
		k =  2 + 2*self.n_breakpoints # No. model parameters
		rss = self.residual_sum_squares
		self.bic = n * np.log(rss/n) + k * np.log(n)

		


class Fit:

	def __init__(self, xx, yy, start_values, max_iterations=30, tolerance=10**-5,
		min_distance_between_breakpoints=0.01, min_distance_to_edge=0.02, verbose=True):

		self.xx = self._validate_list_of_numbers(xx, "xx", min_length=3)
		self.yy = self._validate_list_of_numbers(yy, "yy", min_length=3)

		self.max_iterations = self._validate_integer(max_iterations, "max_iterations")	
		self.tolerance = self._validate_number(tolerance, "tolerance")
		# In terms of proportion of the data range
		self.min_distance_between_breakpoints = self._validate_number(min_distance_between_breakpoints, "min_distance_between_breakpoints")
		# In terms of quantiles
		self.min_distance_to_edge = self._validate_number(min_distance_to_edge, "min_distance_to_edge")
		self.verbose = self._validate_boolean(verbose, "verbose")

		self.start_values = self._validate_start_values(start_values)
		self.n_breakpoints = len(start_values)		

		self.fit_history = []

		self.converged = False
		self.stop = False

		self.davies = None

		if self.verbose:
			print("Input data seems okay")		

		self.fit()


	def _validate_boolean(self, var, var_name):
		if isinstance(var, bool):
			return var
		else:
			raise ValueError("{} must be a Boolean: True or False".format(var_name))


	def _validate_integer(self, var, var_name):
		if isinstance(var, int):
			return var
		else:
			raise ValueError("{} must be an Integer".format(var_name))

	def _validate_number(self, var, var_name):
		if isinstance(var, float) or isinstance(var, int):
			return var
		else:
			raise ValueError("{} must be a Float".format(var_name))

	def _validate_list_of_numbers(self, var, var_name, min_length):
		"""
		Allowed types:
			List of integers of floats
			Numpy array of integers or floats
		"""
		value_error_text = "{} must be a list of numbers with minimum length {}".format(var_name, min_length)
		# If its a list, convert it to a numpy array
		if isinstance(var, list):
			var = np.array(var)

		# If its not a numpy array at this point, raise a value error		
		if not isinstance(var, np.ndarray):
			raise ValueError(value_error_text)

		# Check the array has numebrs in it
		if not np.issubdtype(var.dtype, np.number):
			raise ValueError(value_error_text)

		if len(var) < min_length:
			raise ValueError(value_error_text)

		return var


	def _validate_start_values(self, start_values):

		start_values = self._validate_list_of_numbers(start_values, "start_values", min_length=1)

		if not self._are_breakpoint_values_within_range(start_values):
			raise ValueError("Start values are not within allowed range")
		if not self._are_breakpoint_values_far_apart(start_values):
			raise ValueError("Start values are too close together")
		
		return start_values


	def _are_breakpoint_values_within_range(self, breakpoints):

		min_allowed_bp = np.quantile(self.xx, self.min_distance_to_edge)
		max_allowed_bp = np.quantile(self.xx, 1-self.min_distance_to_edge)

		for bp in breakpoints:
			if bp <= min_allowed_bp or bp >= max_allowed_bp:
				return False				
		return True

	def _are_breakpoint_values_far_apart(self, breakpoints):

		min_distance = np.diff(np.sort(breakpoints))

		# numpy ptp gives the range of the data, closeness realtive to the range
		min_distance_allowed = self.min_distance_between_breakpoints * np.ptp(self.xx)

		if (min_distance <= min_distance_allowed).any():
			return False
		return True


	def fit(self):
		# Run the breakpoint iterative procedure
		if self.verbose:
			print("Running algorithm . . . ")
		while not self.stop:
			# Do the fit


			if len(self.fit_history) == 0:
				current_breakpoints = self.start_values
			else:
				current_breakpoints = self.fit_history[-1]["next_breakpoints"]


			IteratedFit = NextBreakpoints(self.xx, self.yy, current_breakpoints)
			
			next_fit_details = vars(IteratedFit)
			self.fit_history.append(next_fit_details)

			print(next_fit_details["next_breakpoints"])

			self.stop_or_not()

		# Get final davies result
		self.davies = davies.davies_test(self.xx, self.yy)
		
	def stop_or_not(self):

		# Stop if the breakpoints are out of range
		if not self._are_breakpoint_values_within_range(self.fit_history[-1]["next_breakpoints"]):
			if self.verbose:
				print("Breakpoints have diverged out of the data range. Stopping")
			self.stop = True
		
		# Stop if the breakpoints are too close together
		if not self._are_breakpoint_values_far_apart(self.fit_history[-1]["next_breakpoints"]):
			if self.verbose:
				print("Breakpoints have become too close together. Stopping")
			self.stop = True


		# Stop if maximum iterations reached
		if len(self.fit_history)>self.max_iterations:
			if self.verbose:
				print("Max iterations reached. Stopping.")
			self.stop = True

		# Stop if tolerance reached - small change between this and last breakpoints
		if len(self.fit_history) > 1:
			breakpoint_differences = self.fit_history[-2]["next_breakpoints"] - self.fit_history[-1]["next_breakpoints"]
			if np.max(np.abs(breakpoint_differences)) <= self.tolerance:
				if self.verbose:
					print("Algorithm has converged on breakpoint values. Stopping.")
				self.stop = True
				self.converged = True

		# Stop if the algorithm is iterating back to previous values, within tolerance
		if len(self.fit_history) > 2:
			breakpoint_two_step_differences = self.fit_history[-3]["next_breakpoints"] - self.fit_history[-1]["next_breakpoints"]
			if np.max(np.abs(breakpoint_two_step_differences)) <= self.tolerance:
				if self.verbose:
					print("Algorithm is iterating between breakpoint values. Stopping.")
				self.stop = True
				self.converged = True


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
		for bp_i in range(self.n_breakpoints):
			bp_ci = self.estimates["breakpoint{}".format(bp_i+1)]["confidence_interval"]
			plt.axvspan(bp_ci[0], bp_ci[1], alpha=0.1)

	def plot(self, **kwargs):
		"""
		Plot the data, fit, breakpoint poisitons and breakpoint confidence intervals
		"""
		self.plot_data()
		self.plot_fit()
		self.plot_breakpoints()
		self.plot_breakpoint_confidence_intervals()


	def summary(self):
		header = "\n{:^70}\n".format("Breakpoint Regression Results")

		line_length=100
		double_line = "=" * line_length + "\n"
		single_line = "-" * line_length + "\n"

		# Overview
		n_obs = len(self.xx)
		n_model_params = 2 + 2*self.n_breakpoints
		dof = n_obs - n_model_params
		no_obs_text = "{:<20} {:>20}\n".format("No. Observations", n_obs)
		no_model_parameters_text = "{:<20} {:>20}\n".format("No. Model Parameters", n_model_params)
		dof_text = "{:<20} {:>20}\n".format("Degrees of Freedom", dof)
		rss_text = "{:<20} {:>20.6}\n".format("Res. Sum of Squares", self.residual_sum_squares)
		tss_text = "{:<20} {:>20.6}\n".format("Total Sum of Squares", self.total_sum_squares)
		r_2_text = "{:<20} {:>20.6f}\n".format("R Squared", self.r_squared)
		adj_r_2_text = "{:<20} {:>20.6f}\n".format("Adjusted R Squared", self.adjusted_r_squared)

		overview = double_line + no_obs_text + no_model_parameters_text + dof_text + rss_text + tss_text + r_2_text + adj_r_2_text + double_line

		# Table of results

		table_header_template = "{:<15} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12}\n"

		table_header = table_header_template.format("", "Estimate", "Std Err", "t", "P>|t|", "[0.025", "0.975]")
		#print ("{:<8} {:<15} {:<10}".format( name, age, perc))

		table_row_template = "{:<15} {:>12.6} {:>12.3} {:>12.5} {:>12.3} {:>12.5} {:>12.5}\n"

		table_contents = ""

		beta_names = ["beta{}".format(i+1) for i in range(self.n_breakpoints)]
		bp_names = ["breakpoint{}".format(i+1) for i in range(self.n_breakpoints)]

		model_estimator_names = ["const", "alpha1"] + beta_names + bp_names

		for est_name in model_estimator_names:
			estimator_row = table_row_template.format(est_name, self.estimates[est_name]["estimate"], self.estimates[est_name]["se"], 
			self.estimates[est_name]["t_stat"], self.estimates[est_name]["p_t"], 
			self.estimates[est_name]["confidence_interval"][0], self.estimates[est_name]["confidence_interval"][1])
			table_contents += estimator_row

		table_contents += single_line

		table_contents += "These alphas(gradients of segments) are estimated from betas(change in gradient)\n"

		alpha_names = ["alpha{}".format(alpha_i+1) for alpha_i in range(1, self.n_breakpoints+1)]

		table_contents += single_line

		for est_name in alpha_names:
			estimator_row = table_row_template.format(est_name, self.estimates[est_name]["estimate"], self.estimates[est_name]["se"], 
			self.estimates[est_name]["t_stat"], self.estimates[est_name]["p_t"], 
			self.estimates[est_name]["confidence_interval"][0], self.estimates[est_name]["confidence_interval"][1])
			table_contents += estimator_row


		table_contents += double_line

		table = double_line + table_header + single_line + table_contents

		print(header + overview + table)

		print("Alternative hypothesis for breakpoints is not that they equal 0 but that they don't exist - t-stat has no meaning")
		
		if self.davies < 0.05:
			print("Davies test for existence of at least one breakpoint: Null hypothesis of no breakpoint can be ruled out at significance of p<0.05 (p-value is {:.3f})\n".format(self.davies))
		else:
			print("Davies test for existence of at least one breakpoint. Null hypothesis of no breakpoint cannot be ruled out at significance of p<0.05 (p-value is {:.3f})\n".format(self.davies))





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

	print(results.summary())

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

		

	bp_fit = Fit(xx, yy, start_values=[5])


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
	plt.show()



	bp_fit.plot_data()
	bp_fit.plot_fit(color="red", linewidth=4)
	bp_fit.plot_breakpoints()
	bp_fit.plot_breakpoint_confidence_intervals()
	plt.show()



def test_on_data_1c():

	alpha = -4
	beta_1 = -2
	beta_2 = 4
	beta_3 = 1
	intercept = 100
	breakpoint_1 = 7
	breakpoint_2 = 13
	breakpoint_3 = 14

	n_points = 200

	xx = np.linspace(0, 20, n_points)

	yy = intercept + alpha*xx 
	yy += beta_1 * np.maximum(xx - breakpoint_1, 0) 
	yy += beta_2 * np.maximum(xx - breakpoint_2, 0)  
	yy += beta_3 * np.maximum(xx - breakpoint_3, 0)
	yy += np.random.normal(size=n_points)


	bp_fit = Fit(xx, yy, start_values=[5, 10, 16])

	"""
	bp_fit.plot_data()
	bp_fit.plot_fit(color="red", linewidth=4)
	bp_fit.plot_breakpoints()
	bp_fit.plot_breakpoint_confidence_intervals()


	print("The fit data: ", bp_fit.__dict__)


	plt.show()
	plt.close()

	bp_fit.plot_breakpoint_history()
	plt.legend()
	plt.show()

	"""

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
	test_on_data_1c()