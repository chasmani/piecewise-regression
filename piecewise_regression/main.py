
import math

import numpy as np
import statsmodels.api as sm
import scipy.stats
import matplotlib.pyplot as plt

try:
	import piecewise_regression.davies as davies 
	import piecewise_regression.r_squared_calc as r_squared_calc
	from piecewise_regression.data_validation import validate_number, validate_boolean, validate_integer, validate_list_of_numbers
except:
	import davies
	import r_squared_calc
	from data_validation import validate_number, validate_boolean, validate_integer, validate_list_of_numbers


class NextBreakpoints:
	"""
	One iteration of Muggeo's segmented regression algorithm
	Gets the next breakpoints. 
	Also gets interesting statistics etc
	"""

	def __init__(self, xx, yy, current_breakpoints):

		self.xx = validate_list_of_numbers(xx, "xx", min_length=3)
		self.yy = validate_list_of_numbers(yy, "yy", min_length=3)
		self.current_breakpoints = validate_list_of_numbers(current_breakpoints, "current_breakpoints", min_length=1)
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


class Muggeo:
	"""
	Muggeo's iterative segmented regression method
	Raises an error if input variables are the wrong type
	"""


	def __init__(self, xx, yy, start_values, max_iterations=30, tolerance=10**-5,
		min_distance_between_breakpoints=0.01, min_distance_to_edge=0.02, verbose=True):

		self.xx = validate_list_of_numbers(xx, "xx", min_length=3)
		self.yy = validate_list_of_numbers(yy, "yy", min_length=3)

		self.max_iterations = validate_integer(max_iterations, "max_iterations")	
		self.tolerance = validate_number(tolerance, "tolerance")
		# In terms of proportion of the data range
		self.min_distance_between_breakpoints = validate_number(min_distance_between_breakpoints, "min_distance_between_breakpoints")
		# In terms of quantiles
		self.min_distance_to_edge = validate_number(min_distance_to_edge, "min_distance_to_edge")
		self.verbose = validate_boolean(verbose, "verbose")

		self.start_values = self._validate_start_values(start_values)
		self.n_breakpoints = len(start_values)

		self.fit_history = []
		self.best_fit = None

		self.converged = False
		self.stop = False
		# Note records the rason why the algorithm stopped
		self.stop_reason = None

		self.davies = None

		if self.verbose:
			print("Muggeo input data seems okay")		

		self.fit()


	def _validate_start_values(self, start_values):

		start_values = validate_list_of_numbers(start_values, "start_values", min_length=1)

		if not self._are_breakpoint_values_within_range(start_values):
			self.stop = True
			self.stop_reason = "Algorithm was stopped because breakpoint start_values were not within allowed range"
		if not self._are_breakpoint_values_far_apart(start_values):
			self.stop = True
			self.stop_reason = "Algorithm was stopped as breakpoint start_values were too close together"
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
			print("Running Muggeo's iterative algorithm . . . ")
		while not self.stop:
			# Do the fit
			if len(self.fit_history) == 0:
				current_breakpoints = self.start_values
			else:
				current_breakpoints = self.fit_history[-1].next_breakpoints

			IteratedFit = NextBreakpoints(self.xx, self.yy, current_breakpoints)
			
			self.fit_history.append(IteratedFit)
			
			self.stop_or_not()

		# Select the best fit from the fit history by finding the smallest rss
		self.best_fit = min(self.fit_history, key=lambda x:x.residual_sum_squares)


	def stop_or_not(self):

		# Stop if the breakpoints are out of range
		if not self._are_breakpoint_values_within_range(self.fit_history[-1].next_breakpoints):
			self.stop_reason = "Algorithm was stopped as breakpoints diverged out of the data range"
			self.stop = True
		
		# Stop if the breakpoints are too close together
		if not self._are_breakpoint_values_far_apart(self.fit_history[-1].next_breakpoints):
			self.stop_reason = "Algorithm was stopped as breakpoints became too close together"
			self.stop = True

		# Stop if maximum iterations reached
		if len(self.fit_history)>self.max_iterations:
			self.stop_reason = "Algorithm was stopped as maximum iterations reached"
			self.stop = True

		# Stop if tolerance reached - small change between this and last breakpoints
		if len(self.fit_history) > 1:
			breakpoint_differences = self.fit_history[-2].next_breakpoints - self.fit_history[-1].next_breakpoints
			if np.max(np.abs(breakpoint_differences)) <= self.tolerance:
				self.stop_reason = "Algorithm converged on breakpoint values"
				self.stop = True
				self.converged = True

		# Stop if the algorithm is iterating back to previous values, within tolerance
		if len(self.fit_history) > 2:
			breakpoint_two_step_differences = self.fit_history[-3].next_breakpoints - self.fit_history[-1].next_breakpoints
			if np.max(np.abs(breakpoint_two_step_differences)) <= self.tolerance:
				self.stop_reason = "Algorithm converged on breakpoint values"
				self.stop = True
				self.converged = True

	def plot_data(self, **kwargs):
		"""
		Plot the data as a scatter plot
		Passes any kwargs to the matplotlib scatter function, e.g. color="red"
		"""
		plt.scatter(self.xx, self.yy, **kwargs)

	def plot_fit(self, fit_data="best", **kwargs):
		"""
		Plot the fitted model
		Passes any kwargs to the matplotlib plot function, e.g. color="red"
		"""
		if not self.converged:
			print("Algorithm didn't converge. Plotting best fit anyway.")

		# Get the final results from the fitted model variables
		# Params are in terms of [intercept, alpha, betas, gammas]
		
		if fit_data == "best":
			final_params = self.best_fit["raw_params"]
			breakpoints = self.best_fit["next_breakpoints"]
		else:
			final_params = fit_data["raw_params"]
			breakpoints = fit_data["next_breakpoints"]
		
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

	def plot_breakpoints(self, fit_data="best", **kwargs):
		"""
		Plot the breakpoint locations
		Passes kwargs to the matplotlib function, e.g. color="red"
		"""
		if not self.converged:
			print("Algorithm didn't converge. Plotting best fit anyway.")

		if fit_data == "best":
			breakpoints = self.best_fit["next_breakpoints"]			
		else: 
			breakpoints = fit_data["next_breakpoints"]

		for bp in breakpoints:
			plt.axvline(bp, **kwargs)

	def plot_breakpoint_confidence_intervals(self, fit_data="best", **kwargs):
		"""
		Plot the breakpoint cis as shaded regions
		"""
	
		if not self.converged:
			print("Algorithm didn't converge. Plotting breakpoint confidence intervals anyway")
		
		if fit_data == "best":
			estimates = estimates			
		else: 
			estimates = fit_data["estimates"]

		for bp_i in range(self.n_breakpoints):
			bp_ci = estimates["breakpoint{}".format(bp_i+1)]["confidence_interval"]
			plt.axvspan(bp_ci[0], bp_ci[1], alpha=0.1)

	def plot_breakpoint_history(self, **kwargs):
		"""
		Plot the history of the breakpoints as they iterate
		"""
		# Change the format of the breakpoint histories to seperate lists for each breakpoint

		breakpoint_history = [self.start_values]

		for fit_details in self.fit_history:
			breakpoint_history.append(fit_details["next_breakpoints"])

		breakpoint_history_series = zip(*breakpoint_history)
		# Plot a history for each breakpoint 
		count = 0
		for bh in breakpoint_history_series:
			count += 1
			plt.plot(range(0, len(bh)), bh, label="Breakpoint {}".format(count), **kwargs)
			plt.xlabel("Iteration")
			plt.ylabel("Breakpoint")

	def plot_rss_history(self, **kwargs):

		rss_history = []
		for fit_details in self.fit_history:
			rss_history.append(fit_details["residual_sum_squares"])

		plt.plot(range(1, len(rss_history) +1), rss_history, **kwargs)
		plt.xlabel("Iteration")
		plt.ylabel("Residual Sum of Squares")


	def plot(self, **kwargs):
		"""
		Plot the data, fit, breakpoint poisitons and breakpoint confidence intervals
		"""
		self.plot_data()
		self.plot_fit()
		self.plot_breakpoints()
		self.plot_breakpoint_confidence_intervals()

	def summary(self):
		if self.converged:
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
			rss_text = "{:<20} {:>20.6}\n".format("Res. Sum of Squares", self.best_fit.residual_sum_squares)
			tss_text = "{:<20} {:>20.6}\n".format("Total Sum of Squares", self.best_fit.total_sum_squares)
			r_2_text = "{:<20} {:>20.6f}\n".format("R Squared", self.best_fit.r_squared)
			adj_r_2_text = "{:<20} {:>20.6f}\n".format("Adjusted R Squared", self.best_fit.adjusted_r_squared)

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

			estimates = self.best_fit.estimates

			for est_name in model_estimator_names:
				estimator_row = table_row_template.format(est_name, estimates[est_name]["estimate"], estimates[est_name]["se"], 
				estimates[est_name]["t_stat"], estimates[est_name]["p_t"], 
				estimates[est_name]["confidence_interval"][0], estimates[est_name]["confidence_interval"][1])
				table_contents += estimator_row

			table_contents += single_line

			table_contents += "These alphas(gradients of segments) are estimated from betas(change in gradient)\n"

			alpha_names = ["alpha{}".format(alpha_i+1) for alpha_i in range(1, self.n_breakpoints+1)]

			table_contents += single_line

			for est_name in alpha_names:
				estimator_row = table_row_template.format(est_name, estimates[est_name]["estimate"], estimates[est_name]["se"], 
				estimates[est_name]["t_stat"], estimates[est_name]["p_t"], 
				estimates[est_name]["confidence_interval"][0], estimates[est_name]["confidence_interval"][1])
				table_contents += estimator_row


			table_contents += double_line

			table = double_line + table_header + single_line + table_contents

			print(header + overview + table)

			print("Alternative hypothesis for breakpoints is not that they equal 0 but that they don't exist - t-stat has no meaning")

			"""
			
			if self.davies < 0.05:
				print("Davies test for existence of at least one breakpoint: Null hypothesis of no breakpoint can be ruled out at significance of p<0.05 (p-value is {:.3f})\n".format(self.davies))
			else:
				print("Davies test for existence of at least one breakpoint. Null hypothesis of no breakpoint cannot be ruled out at significance of p<0.05 (p-value is {:.3f})\n".format(self.davies))
			"""

class BootstrapRestarting:

	
	def __init__(self, xx, yy, n_breakpoints, start_values=None, max_iterations=30, tolerance=10**-5,
		min_distance_between_breakpoints=0.01, min_distance_to_edge=0.02, verbose=True,
		n_boot=10):

		self.xx = validate_list_of_numbers(xx, "xx", min_length=3)
		self.yy = validate_list_of_numbers(yy, "yy", min_length=3)

		self.max_iterations = validate_integer(max_iterations, "max_iterations")	
		self.tolerance = validate_number(tolerance, "tolerance")
		# In terms of proportion of the data range
		self.min_distance_between_breakpoints = validate_number(min_distance_between_breakpoints, "min_distance_between_breakpoints")
		# In terms of quantiles
		self.min_distance_to_edge = validate_number(min_distance_to_edge, "min_distance_to_edge")
		self.verbose = validate_boolean(verbose, "verbose")

		self.n_breakpoints = validate_integer(n_breakpoints, "n_breakpoints")

		self.n_boot = validate_integer(n_boot, "n_boot")

		self.start_values = start_values

		self.bootstrap_history = []
		
		self.best_muggeo = None

		self.stop = False

		if self.verbose:
			print("Input data seems okay")		

		self.bootstrap_restarting()

	def bootstrap_restarting(self):

		if self.verbose:
			print("Running bootstrap restarting . . . ")

		# 1. Run it without bootstrapping first of all
		if not self.start_values:
			if self.verbose:
				print("No start_values given, randomly generating . . . ")
			min_allowed_bp = np.quantile(self.xx, self.min_distance_to_edge)
			max_allowed_bp = np.quantile(self.xx, 1-self.min_distance_to_edge)
			start_bps = np.random.uniform(low=min_allowed_bp, high=max_allowed_bp, size=self.n_breakpoints)

		muggeo_fit = Muggeo(self.xx, self.yy, start_bps)
		self.bootstrap_history.append(muggeo_fit)
		if muggeo_fit.converged:
			self.best_muggeo = muggeo_fit

		# Iterate bootstraps
		count = 1
		while not self.stop:
			
			# Best breakpoints are either from best muggeo so far, start values or random
			if self.best_muggeo:
				best_bps = self.best_muggeo.best_fit.next_breakpoints
			elif self.start_values:
				# This might be None, in which case the Muggeo algo will randomly choose some
				best_bps = self.start_values
			else:
				min_allowed_bp = np.quantile(self.xx, self.min_distance_to_edge)
				max_allowed_bp = np.quantile(self.xx, 1-self.min_distance_to_edge)
				best_bps = np.random.uniform(low=min_allowed_bp, high=max_allowed_bp, size=self.n_breakpoints)

				
			# Get some new breakpoint values from a bootstrapped fit
			xx_boot, yy_boot = self.bootstrap_data(self.xx, self.yy)
			bootstrap_fit = Muggeo(xx_boot, yy_boot, best_bps)
			bootstrap_bps = bootstrap_fit.best_fit.next_breakpoints

			# Do a new fit with the new breakpoint values 
			next_muggeo = Muggeo(self.xx, self.yy, bootstrap_bps)
			self.bootstrap_history.append(next_muggeo)

			# If we get a converged answer, see if this new fit is the best or not
			if next_muggeo.converged:
				# If there is already a best_muggeo, see if this one is better
				if self.best_muggeo:
					if next_muggeo.best_fit.residual_sum_squares < self.best_muggeo.best_fit.residual_sum_squares:
						self.best_muggeo = next_muggeo
						print("New best fit with rss ", next_muggeo.best_fit.residual_sum_squares)
				else:
					self.best_muggeo = next_muggeo

			self.stop_or_not_bootstrap()


	def bootstrap_data(self, xx, yy):
		"""
		Non parametric bootstrap, randomly sample data points with replacement
		Return bootstrapped data of same length as oriignal data
		"""

		n = len(xx)
		# Get bootstrap samples as array index locations
		boot_indices = np.random.choice(n, size=n, replace=True)

		xx_boot = xx[boot_indices]
		yy_boot = yy[boot_indices]
		return xx_boot, yy_boot

	def stop_or_not_bootstrap(self):

		# Stop if maximum iterations reached
		if len(self.bootstrap_history)>=self.n_boot:
			if self.verbose:
				print("Max bootstrap iterations reached. Stopping.")
			self.stop = True


class ModelSection:

	def __init__(self, xx, yy, max_breakpoints, n_boot=10, start_values=None, max_iterations=30, tolerance=10**-5,
		min_distance_between_breakpoints=0.01, min_distance_to_edge=0.02, verbose=True):

		self.xx = validate_list_of_numbers(xx, "xx", min_length=3)
		self.yy = validate_list_of_numbers(yy, "yy", min_length=3)

		self.max_iterations = validate_integer(max_iterations, "max_iterations")	
		self.tolerance = validate_number(tolerance, "tolerance")
		# In terms of proportion of the data range
		self.min_distance_between_breakpoints = validate_number(min_distance_between_breakpoints, "min_distance_between_breakpoints")
		# In terms of quantiles
		self.min_distance_to_edge = validate_number(min_distance_to_edge, "min_distance_to_edge")
		self.verbose = validate_boolean(verbose, "verbose")

		self.n_boot = validate_integer(n_boot, "n_boot")

		self.max_breakpoints = validate_integer(max_breakpoints, "max_breakpoints")

		self.models = []
		
		self.stop = False

		self.model_selection()

	def model_selection(self):

		for k in range(1, self.max_breakpoints+1):
			bootstrapped_fit = BootstrapRestarting(self.xx, self.yy, n_breakpoints=k)
			self.models.append(bootstrapped_fit)
			
			if not bootstrapped_fit.best_muggeo:
				break

		for model in self.models:
			print(model.best_muggeo.next_breakpoints, model.best_muggeo.best_fit.bic, model.best_muggeo.best_fit.residual_sum_squares)


	

class Fit:
	"""
	Set up the fit etc

	"""
	pass


	#self.davies = None


if __name__=="__main__":
	alpha = -4
	beta_1 = -2
	intercept = 100
	breakpoint_1 = 7

	n_points = 200

	xx = np.linspace(0, 20, n_points)
	yy = intercept + alpha*xx + beta_1 * np.maximum(xx - breakpoint_1, 0) + np.random.normal(size=n_points)

	
	fit = ModelSection(xx, yy, max_breakpoints=10)
