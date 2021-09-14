

import math

import numpy as np
import statsmodels.api as sm
import scipy.stats
import matplotlib.pyplot as plt

try:
	import piecewise_regression.r_squared_calc as r_squared_calc
	from piecewise_regression.data_validation import validate_positive_number, validate_boolean, validate_positive_integer, validate_list_of_numbers, validate_non_negative_integer
	from piecewise_regression.main import Fit
except:
	import r_squared_calc
	from data_validation import validate_positive_number, validate_boolean, validate_positive_integer, validate_list_of_numbers, validate_non_negative_integer
	from main import Fit

class ModelSelection:
	"""
	Experimental - uses simple BIC based on simple linear model.
	"""

	def __init__(self, xx, yy, max_breakpoints=10, n_boot=20, max_iterations=30, tolerance=10**-5,
		min_distance_between_breakpoints=0.01, min_distance_to_edge=0.02, verbose=True):


		self.models = []
		self.full_models_data = []
		
		self.stop = False

		if verbose==True:
			print("Running fit with n_breakpoint = 0 . . ")

		self.no_breakpoint_fit(xx, yy)

		for k in range(1, max_breakpoints+1):
			if verbose==True:
				print("Running fit with n_breakpoint = {} . . ".format(k))
			bootstrapped_fit = Fit(xx, yy, n_breakpoints=k, verbose=False, 
				n_boot=n_boot, max_iterations=max_iterations, tolerance=tolerance,
				min_distance_between_breakpoints=min_distance_between_breakpoints, min_distance_to_edge=min_distance_to_edge)
			fit_summary = bootstrapped_fit.get_results()
			fit_summary["n_breakpoints"] = k
			self.models.append(fit_summary)
			self.full_models_data.append(bootstrapped_fit)

		self.summary()


	def summary(self):
		
		header = "\n{:^70}\n".format("Breakpoint Model Comparision Results")

		line_length=100
		double_line = "=" * line_length + "\n"
		single_line = "-" * line_length + "\n"

		table_header_template = "{:<15} {:>12} {:>12} {:>12} \n"
		table_header = table_header_template.format("n_breakpoints", "BIC", "converged", "RSS")		
		table_row_template = "{:<15} {:>12.5} {:>12} {:>12.5} \n"

		table_contents = header
		table_contents += double_line

		table_contents += table_header
		table_contents += single_line

		for model in self.models:

			if model["converged"]:
				model_row = table_row_template.format(model["n_breakpoints"], model["bic"], str(model["converged"]), model["rss"])
			else:
				model_row = table_row_template.format(model["n_breakpoints"], "", str(model["converged"]), "")

			table_contents += model_row

		print(table_contents)

		print("Minimum BIC (Bayesian Information Criterion) suggests the best model")


	def no_breakpoint_fit(self, xx, yy):

		Z = np.array([xx])
		Z = Z.T
		Z = sm.add_constant(Z, has_constant='add')
		# Basic OLS fit
		results = sm.OLS(endog=np.array(yy), exog=Z).fit()

		# get the predicted values 
		ff = [(results.params[0] + results.params[1] * x) for x in xx]

		# Get Rss
		rss, tss, r_2, adjusted_r_2 = r_squared_calc.get_r_squared(yy, ff, n_params=2)

		# Calcualte BIC
		n = len(xx) # No. data points
		k =  2 # No. parameters
		bic = n * np.log(rss/n) + k * np.log(n)

		fit_data = {
			"bic":bic,
			"n_breakpoints":0,
			"estimates":{},
			"converged":True,
			"rss":rss
		}

		fit_data["estimates"]["const"] = results.params[0]
		fit_data["estimates"]["alpha1"] = results.params[1]

		self.models.append(fit_data)
		self.full_models_data.append(fit_data)


if __name__=="__main__":

	#np.random.seed(2)
	pass