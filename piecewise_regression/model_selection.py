


class ModelSection:
	"""
	Experimental - uses simple BIC based on simple linear model.
	"""

	def __init__(self, xx, yy, n_breakpoints, max_breakpoints=10, n_boot=10, start_values=None, max_iterations=30, tolerance=10**-5,
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
			
			#if not bootstrapped_fit.best_muggeo:
			#	break

		for model in self.models:
			print(model.best_muggeo)

			#print(model.best_muggeo.next_breakpoints, model.best_muggeo.best_fit.bic, model.best_muggeo.best_fit.residual_sum_squares)

