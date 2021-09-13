

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import math

import os, sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from piecewise_regression import Fit


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

	bp_fit.plot_best_muggeo_breakpoint_history()
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




	

if __name__=="__main__":

	np.random.seed(0)
	test_on_data_1b()