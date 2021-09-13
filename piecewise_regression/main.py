
import math

import numpy as np
import statsmodels.api as sm
import scipy.stats
import matplotlib.pyplot as plt

try:
	import piecewise_regression.davies as davies 
	import piecewise_regression.r_squared_calc as r_squared_calc
	from piecewise_regression.data_validation import validate_positive_number, validate_boolean, validate_positive_integer, validate_list_of_numbers
except:
	import davies
	import r_squared_calc
	from data_validation import validate_positive_number, validate_boolean, validate_positive_integer, validate_list_of_numbers


class NextBreakpoints:
	"""
	One iteration of Muggeo's segmented regression algorithm
	Gets the next breakpoints. 
	Also gets interesting statistics etc
	"""
	def __init__(self, xx, yy, current_breakpoints):

		# Data validation done at a higher level
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

		if self.verbose:
			print("Iterating fit . . . ")
			print("Current breakpoints are: ", self.current_breakpoints)
			print("Next breakpoints are: ", self.next_breakpoints)


	def calculate_all_estimates(self):
		"""
		Save all params in self.estimates dictionary
		"""
		params = self.raw_params
	
		# Extract the exstimates form the correct locations in the params returned from statsmodels
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
		"""
		Get the standard errors for the alphas (gradients of segments)
		"""
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
		"""
		Get the standard errors of the breakpoints
		bp = gamma/beta + bp_0
		Variance of bp estimator found using ratio/delta method
		See e.g. Muggeo (2003) for clarification
		"""
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

			# From Muggeo (2003). The sign before the covariance term is opposite to Muggeo, this is because the gamma is defined with the opposite sign
			# The calculation is equivalent to Muggeos.
			bp_var = (gamma_var + beta_var * (gamma/beta)**2 - 2*(gamma/beta)*gamma_beta_covar)/(beta**2)
			bp_vars.append(bp_var)

		bp_ses = np.sqrt(bp_vars)
		return bp_ses


	def get_const_standard_error(self):
		# Covariance matrix is [c, alpha, betas, gammas]
		# Constant variance is just the top left cell in covariance matrix
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
			if "breakpoint" in estimator_name or "beta" in estimator_name:
				details["t_stat"] = "-"
				details["p_t"] = "-"
			else:
				t_stat = details["estimate"]/details["se"]
				p_t = scipy.stats.t.sf(np.abs(t_stat), dof)*2
				details["t_stat"] = t_stat
				details["p_t"] = p_t		


	def get_predicted_yy(self):
		"""
		Get the model predictions for each of the xx data points
		"""
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
		"""
		Calculate R squared from fitted model
		Uses an imported function from another file
		"""
		yy_predicted = self.get_predicted_yy()
		n_params = 2 * self.n_breakpoints + 2

		self.residual_sum_squares, self.total_sum_squares, self.r_squared, self.adjusted_r_squared = r_squared_calc.get_r_squared(self.yy, yy_predicted, n_params)

	def calculate_bayesian_information_criterion(self):
		"""
		Assuming normal noise, uses the standard version for OLS models. 
		I beleive this holds for brekpoint regression models, because the BIC is based on the likelihood of the data
		given the model. That likelihood function won't include the breakpoint values - it just depends on distances of the data to the fitted model predictions
		Also depends on the error in the noise term, should work as long as the noise is constant.  
		"""
		n = len(self.xx) # No. data points
		k =  2 + 2*self.n_breakpoints # No. model parameters
		rss = self.residual_sum_squares
		self.bic = n * np.log(rss/n) + k * np.log(n)


class Muggeo:
	"""
	Muggeo's iterative segmented regression method
	Simple version. Errors are handled at a higher level in the Fit object
	See Muggeo (2003) for more information
	"""
	def __init__(self, xx, yy, start_values, max_iterations=30, tolerance=10**-5, verbose=True):

		if self.verbose:
			print("Instantiating Muggeo . . . with start_values ={}".format(start_values))

		# validation is done at a higher level
		self.xx = xx
		self.yy = yy
		self.start_values = start_values
		self.max_iterations = max_iterations	
		self.tolerance = tolerance
		self.verbose = verbose


		self.stop = False

		self.n_breakpoints = len(start_values)
		self.fit_history = []
		self.best_fit = None
		self.converged = False
		
		# Records the reason why the algorithm stopped
		self.stop_reason = None

		self.fit()


	def fit(self):
		# Run the breakpoint iterative procedure
		if self.verbose:
			print("Running Muggeo's iterative algorithm . . . ")
		while not self.stop:
			# Do the fit
			# Get the current breakpoints. If first fit then use start_vakues, otherwise use fit history
			if len(self.fit_history) == 0:
				current_breakpoints = self.start_values
			else:
				current_breakpoints = self.fit_history[-1].next_breakpoints
			try:
				IteratedFit = NextBreakpoints(self.xx, self.yy, current_breakpoints)
				self.fit_history.append(IteratedFit)
				self.stop_or_not()
			except Exception as e:
				self.stop=True
				self.stop_reason = "Error encountered: " + str(e)

		# Select the best fit from the fit history by finding the smallest rss
		if len(self.fit_history) > 0:
			self.best_fit = min(self.fit_history, key=lambda x:x.residual_sum_squares)

	def stop_or_not(self):
		"""
		Stop if it's converegd or max_iterations reached
		"""
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


class Fit:
	"""
	Fit a segmented regression model to data
	Uses bootstrap restarting to avoid local minima
	Requires either n_breakpoints of start_values
	if no start_vaues are given, they are instead uniformly randomly generated across range of data
	Also variabels to control how the fit is run
	"""
	def __init__(self, xx, yy,
		n_breakpoints=None, start_values=None, n_boot=10, verbose=True, 
		max_iterations=30, tolerance=10**-5,
		min_distance_between_breakpoints=0.01, min_distance_to_edge=0.02,
		):

		# Validate all input data
		self.xx = validate_list_of_numbers(xx, "xx", min_length=3)
		self.yy = validate_list_of_numbers(yy, "yy", min_length=3)
		self.n_boot = validate_positive_integer(n_boot, "n_boot")
		self.verbose = validate_boolean(verbose, "verbose")
		self.max_iterations = validate_positive_integer(max_iterations, "max_iterations")
		self.tolerance = validate_positive_number(tolerance, "tolerance")
		self.min_distance_between_breakpoints = validate_positive_number(min_distance_between_breakpoints, "min_distance_between_breakpoints")
		self.min_distance_to_edge = validate_positive_number(min_distance_to_edge, "min_distance_to_edge")

		# We need either start_values or n_breakpoints
		# start_values takes precedence if n_breakpoints doesn't match
		if start_values:
			self.start_values = _validate_start_values(start_values, "start_values")
			self.n_breakpoints = len(self.start_values)
		else:
			if self.n_breakpoints:
				self.n_breakpoints = validate_positive_integer(n_breakpoints, "n_breakpoints")
			else:
				raise ValueError("The Fit algorithm requires either the start_values or n_breakpoints")				

		self.bootstrap_history = []
		self.best_muggeo = None
		self.stop = False

		self.bootstrap_restarting()

		self.davies = davies.davies_test(self.xx, self.yy)


	def _validate_start_values(self, start_values):


		start_values = validate_list_of_numbers(start_values, "start_values", min_length=1)

		if not self._are_breakpoint_values_within_range(start_values):
			self.stop = True
			raise ValueError("start_values are not within allowed range. Try changing min_distance_to_edge")
			
		if not self._are_breakpoint_values_far_apart(start_values):
			self.stop = True
			raise ValueError("start_values are too close together. Try changing min_distance_between_breakpoints")
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

	def _generate_valid_breakpoints(self):

		# Get breakpoints within allowed range
		min_allowed_bp = np.quantile(self.xx, self.min_distance_to_edge)
		max_allowed_bp = np.quantile(self.xx, 1-self.min_distance_to_edge)
		start_values = np.random.uniform(low=min(self.xx), high=max(self.xx), size=self.n_breakpoints)

		# Get breakpoints not too close together
		gen_count = 0
		while not self._are_breakpoint_values_far_apart(start_values) and not self._are_breakpoint_values_within_range(start_values):
			start_values = np.random.uniform(low=min(self.xx), high=max(self.xx), size=self.n_breakpoints)
			gen_count += 1
			# Stop the genreation after 50 attempts. Almost certainly this will be because the disatnce between breakpoints is too small
			if gen_count == 100:
				raise ValueError("Unable to generate random breakpoints that are far enough apart. Change min_distance_between_breakpoints")
		return start_values

	def bootstrap_restarting(self):
		
		# Do a first fit with the start_values given. Use random ones if no start_values
		if self.start_values:
			start_bps = self.start_values
		else:
			start_bps = self._generate_valid_breakpoints()

		muggeo_fit = Muggeo(xx=self.xx, yy=self.yy, start_values=start_bps, 
			max_iterations=self.max_iterations, tolerance = self.tolerance, verbose=self.verbose)

		self.bootstrap_history.append(muggeo_fit)

		# Start off with the "best_muggeo" even if it didn't converge
		self.best_muggeo = muggeo_fit

		# Iterate bootstraps
		for i in range(self.n_boot):

			print("Looping restarting . . ")
			# Best breakpoints are either from best converged muggeo so far, or start values, or randomly generated
			if self.best_muggeo.converged:
				best_bps = self.best_muggeo.best_fit.next_breakpoints
			else:
				if self.start_values:
					best_bps = self.start_values
				else:
					best_bps = self._generate_valid_breakpoints()
				
			# Get some new breakpoint values from a bootstrapped fit
			# Non parametric bootstrap by resampling from data
			xx_boot, yy_boot = self.bootstrap_data(self.xx, self.yy)
			bootstrap_fit = Muggeo(xx_boot, yy_boot, start_values=best_bps, 
				max_iterations=self.max_iterations, tolerance = self.tolerance, verbose=self.verbose)
			if bootstrap_fit.converged:
				bootstrap_bps = bootstrap_fit.best_fit.next_breakpoints
			else:
				# Give it something - even though these breakpoints are already run
				# prefer this to using breaks in the for loop
				bootstrap_bps = best_bps

			# Do a new fit with the new breakpoint values 
			next_muggeo = Muggeo(self.xx, self.yy, start_values=bootstrap_bps, 
				max_iterations=self.max_iterations, tolerance = self.tolerance, verbose=self.verbose)
			self.bootstrap_history.append(next_muggeo)

			# If we get a converged answer, see if this new fit is the best or not
			if next_muggeo.converged:
				print("Converged")
				# If there is already a converged best_muggeo, see if this one is better
				if self.best_muggeo.converged:
					if next_muggeo.best_fit.residual_sum_squares < self.best_muggeo.best_fit.residual_sum_squares:
						self.best_muggeo = next_muggeo
						print("New best fit with rss ", next_muggeo.best_fit.residual_sum_squares)
				# If there is not already a converged best_muggeo, use this fit instead
				else:
					self.best_muggeo = next_muggeo


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
		if not self.best_muggeo.converged:
			print("Algorithm didn't converge. Plotting best fit anyway (not very meaningful)")

		# Get the final results from the fitted model variables
		# Params are in terms of [intercept, alpha, betas, gammas]
		final_params = self.best_muggeo.best_fit.raw_params
		breakpoints = self.best_muggeo.best_fit.next_breakpoints
		
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
		if not self.best_muggeo.converged:
			print("Algorithm didn't converge. No meaningful breakpoint estimates to plot")
		else:
			breakpoints = self.best_muggeo.best_fit.next_breakpoints			
			
			for bp in breakpoints:
				plt.axvline(bp, **kwargs)

	def plot_breakpoint_confidence_intervals(self, **kwargs):
		"""
		Plot the breakpoint cis as shaded regions
		"""
	
		if not self.best_muggeo.converged:
			print("Algorithm didn't converge. No meaningful breakpoint estimates to plot")
		else:			
			estimates = self.best_muggeo.best_fit.estimates			
			
			for bp_i in range(self.best_muggeo.n_breakpoints):
				bp_ci = estimates["breakpoint{}".format(bp_i+1)]["confidence_interval"]
				plt.axvspan(bp_ci[0], bp_ci[1], alpha=0.1)

	def plot_breakpoint_history(self, **kwargs):
		"""
		Plot the history of the breakpoints as they iterate. 
		History of the best muggeo fit.
		"""
		if not self.best_muggeo.converged:
			print("Algorithm didn't converge.")

		# Get the data from the best_muggeo in a form for plotting
		breakpoint_history = [self.best_muggeo.start_values]
		for fit_details in self.best_muggeo.fit_history:
			breakpoint_history.append(fit_details.next_breakpoints)
		breakpoint_history_series = zip(*breakpoint_history)
		
		# Plot a history for each breakpoint 
		count = 0
		for bh in breakpoint_history_series:
			count += 1
			plt.plot(range(0, len(bh)), bh, label="Breakpoint {}".format(count), **kwargs)
			plt.xlabel("Iteration")
			plt.ylabel("Breakpoint")

	def plot(self):
		"""
		Plot the data, fit, breakpoint positions and breakpoint confidence intervals
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
		rss_text = "{:<20} {:>20.6}\n".format("Res. Sum of Squares", self.best_muggeo.best_fit.residual_sum_squares)
		tss_text = "{:<20} {:>20.6}\n".format("Total Sum of Squares", self.best_muggeo.best_fit.total_sum_squares)
		r_2_text = "{:<20} {:>20.6f}\n".format("R Squared", self.best_muggeo.best_fit.r_squared)
		adj_r_2_text = "{:<20} {:>20.6f}\n".format("Adjusted R Squared", self.best_muggeo.best_fit.adjusted_r_squared)
		converged_text = "{:<20} {:>20s}".format("Converged: ", str(self.best_muggeo.converged))
		if not self.best_muggeo.converged:
			converged_text += "  (WARNING: NO CONVERGENCE MEANS UNRELIABLE ESTIMATES IN THIS SUMMARY)\n"
		else:
			converged_text += "\n"

		overview = double_line + no_obs_text + no_model_parameters_text + dof_text + rss_text + tss_text + r_2_text + adj_r_2_text + converged_text + double_line

		# Table of results

		table_header_template = "{:<15} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12}\n"

		table_header = table_header_template.format("", "Estimate", "Std Err", "t", "P>|t|", "[0.025", "0.975]")
		#print ("{:<8} {:<15} {:<10}".format( name, age, perc))

		table_row_template = "{:<15} {:>12.6} {:>12.3} {:>12.5} {:>12.3} {:>12.5} {:>12.5}\n"

		table_contents = ""

		beta_names = ["beta{}".format(i+1) for i in range(self.n_breakpoints)]
		bp_names = ["breakpoint{}".format(i+1) for i in range(self.n_breakpoints)]

		model_estimator_names = ["const", "alpha1"] + beta_names + bp_names

		estimates = self.best_muggeo.best_fit.estimates

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

		print("Davies test for existence of at least 1 breakpoint: p={:.6} (e.g. p<0.05 means reject null hypothesis of no breakpoints at 5% significance)".format(self.davies))
		print("\n\n")


if __name__=="__main__":

	#np.random.seed(2)
	np.random.seed(2)
	alpha = -4
	beta_1 = -2
	intercept = 100
	breakpoint_1 = 17

	n_points = 200

	xx = np.linspace(10, 30, n_points)
	#yy = intercept + alpha*xx + beta_1 * np.maximum(xx - breakpoint_1, 0) + np.random.normal(size=n_points)
	yy = intercept + alpha*xx + beta_1
	
	print(xx, yy)


	xx_str = [str(x) for x in xx]

	print(", ". join(xx_str))

	yy_str = [str(y) for y in yy]

	print(", ". join(yy_str))


	

	fit = Fit(xx, yy, n_breakpoints=1)


	fit.summary()

	fit.plot()
	plt.show()

	