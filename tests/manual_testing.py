

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import math

import os, sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from piecewise_regression import Fit
from piecewise_regression import ModelSelection


def on_data_1():

	alpha = -4
	beta_1 = -2
	intercept = 100
	breakpoint_1 = 7

	n_points = 200

	xx = np.linspace(0, 20, n_points)
	yy = intercept + alpha*xx + beta_1 * np.maximum(xx - breakpoint_1, 0) + np.random.normal(size=n_points)

		

	pw_fit = Fit(xx, yy, start_values=[5])


	print("p-value is ", pw_fit.davies)

	pw_results = pw_fit.get_results()
	pw_estimates = pw_results["estimates"]
	print(pw_results)

	print(pw_estimates)

	pw_bootstrap_history = pw_fit.bootstrap_history
	print(pw_bootstrap_history)



	#print(bp_fit.breakpoint_history)

	#bp_fit.plot_data()
	#plt.show()


def on_data_1b():

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

	bp_fit.plot_best_muggeo_breakpoint_history()
	plt.show()

	bp_fit.plot_data()
	bp_fit.plot_fit(color="red", linewidth=4)
	bp_fit.plot_breakpoints()
	bp_fit.plot_breakpoint_confidence_intervals()
	plt.show()



def on_data_1c():

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

	bp_fit.summary()


	
	bp_fit.plot_data()
	bp_fit.plot_fit(color="red", linewidth=4)
	bp_fit.plot_breakpoints()
	bp_fit.plot_breakpoint_confidence_intervals()


	print("The fit data: ", bp_fit.__dict__)


	plt.show()
	plt.close()

	bp_fit.plot_best_muggeo_breakpoint_history()
	plt.legend()
	plt.show()
	plt.close()

	bp_fit.plot_bootstrap_restarting_history()
	plt.legend()
	plt.show()
	plt.close()



def model_selection_1():

	alpha = -4
	beta_1 = -2
	intercept = 100
	breakpoint_1 = 17

	n_points = 100

	xx = np.linspace(10, 30, n_points)
	yy = intercept + alpha*xx + beta_1 * np.maximum(xx - breakpoint_1, 0) + np.random.normal(size=n_points)
	#yy = intercept + alpha*xx + beta_1 + np.random.normal(size=n_points)
	

	xx_str = [str(x) for x in xx]
	yy_str = [str(y) for y in yy]	

	ms = ModelSelection(xx, yy, max_breakpoints=6)


def model_selection_2():


	alpha = -4
	beta_1 = -4
	beta_2 = 4
	intercept = 100
	breakpoint_1 = 7
	breakpoint_2 = 12

	n_points = 200

	xx = np.linspace(0, 20, n_points)

	yy = intercept + alpha*xx + beta_1 * np.maximum(xx - breakpoint_1, 0) + beta_2 * np.maximum(xx-breakpoint_2, 0)  + np.random.normal(size=n_points)


	ms = ModelSelection(xx, yy)


def fit_3_check_this_makes_sense():

	np.random.seed(0)

	alpha = 10
	beta_1 = -8
	beta_2 = -6
	beta_3 = 10
	intercept = 100
	breakpoint_1 = 7
	breakpoint_2 = 10
	breakpoint_3 = 14

	n_points = 200

	xx = np.linspace(0, 20, n_points)

	yy = intercept + alpha*xx 
	yy += beta_1 * np.maximum(xx - breakpoint_1, 0) 
	yy += beta_2 * np.maximum(xx - breakpoint_2, 0)  
	yy += beta_3 * np.maximum(xx - breakpoint_3, 0)
	yy += np.random.normal(size=n_points)


	pr = Fit(xx, yy, n_breakpoints=2)
	pr.plot()
	plt.show()

	pr3 = Fit(xx, yy, n_breakpoints=3)
	pr3.plot()
	plt.show()

	pr4 = Fit(xx, yy, n_breakpoints=4)
	pr4.plot()
	plt.show()

	ms = ModelSelection(xx, yy, max_breakpoints=6)


def fit_with_initally_diverging():
	np.random.seed(2)

	alpha = 10
	beta_1 = -8
	beta_2 = 3
	beta_3 = 10
	intercept = 100
	breakpoint_1 = 7
	breakpoint_2 = 10
	breakpoint_3 = 14

	n_points = 200

	xx = np.linspace(0, 20, n_points)

	yy = intercept + alpha*xx 
	yy += beta_1 * np.maximum(xx - breakpoint_1, 0) 
	yy += beta_2 * np.maximum(xx - breakpoint_2, 0)  
	yy += beta_3 * np.maximum(xx - breakpoint_3, 0)
	yy += np.random.normal(size=n_points)


	pr = Fit(xx, yy, n_breakpoints=2)
	print(pr.summary)

def fit_with_initially_diverging_start_values():

	np.random.seed(0)

	alpha = 10
	beta_1 = -8
	beta_2 = 3
	beta_3 = 10
	intercept = 100
	breakpoint_1 = 7
	breakpoint_2 = 10
	breakpoint_3 = 14

	n_points = 200

	xx = np.linspace(0, 20, n_points)

	yy = intercept + alpha*xx 
	yy += beta_1 * np.maximum(xx - breakpoint_1, 0) 
	yy += beta_2 * np.maximum(xx - breakpoint_2, 0)  
	yy += beta_3 * np.maximum(xx - breakpoint_3, 0)
	yy += np.random.normal(size=n_points)


	pr = Fit(xx, yy, start_values=[2.15646833, 0.98300926], n_boot=20)
	pr.summary()

def fit_with_initially_diverging_start_values_b():

	np.random.seed(0)

	alpha = 10
	beta_1 = -8
	beta_2 = 3
	beta_3 = 10
	intercept = 100
	breakpoint_1 = 7
	breakpoint_2 = 10
	breakpoint_3 = 14

	n_points = 200

	xx = np.linspace(0, 20, n_points)

	yy = intercept + alpha*xx 
	yy += beta_1 * np.maximum(xx - breakpoint_1, 0) 
	yy += beta_2 * np.maximum(xx - breakpoint_2, 0)  
	yy += beta_3 * np.maximum(xx - breakpoint_3, 0)
	yy += np.random.normal(size=n_points)


	pr = Fit(xx, yy, start_values=[1.2, 0.53], n_boot=25)
	pr.summary()

def fit_with_straight_line():

	np.random.seed(0)

	alpha = 10
	intercept = 100

	n_points = 200

	xx = np.linspace(0, 20, n_points)

	yy = intercept + alpha*xx 
	yy += np.random.normal(size=n_points)


	pr = Fit(xx, yy, n_breakpoints=0, n_boot=25)
	pr.summary()


def model_comparision_straight_line():

	np.random.seed(0)

	alpha = 10
	intercept = 100

	n_points = 200

	xx = np.linspace(0, 20, n_points)

	yy = intercept + alpha*xx 
	yy += np.random.normal(size=n_points)


	ms = ModelSelection(xx, yy, max_breakpoints=6)


if __name__=="__main__":

	on_data_1()