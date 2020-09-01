
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

def generate_breakpoint_data():

	alpha = -4
	beta_1 = -2
	intercept = 100
	breakpoint_1 = 7

	n_points = 200

	xx = np.linspace(0, 20, n_points)

	yy = intercept + alpha*xx + beta_1 * np.maximum(xx - breakpoint_1, 0) + np.random.normal(size=n_points)

	return xx, yy

	"""
	plt.scatter(xx, yy, s=2)
	plt.show()
	"""


def breakpoint_fit(zz, yy, current_breakpoint_1=10):

	A = zz
	B = (zz - current_breakpoint_1) * np.heaviside(zz- current_breakpoint_1, 1) 
	C = np.heaviside(zz - current_breakpoint_1, 1)
	
	Z = np.array([A, B , C]).T

	Z = sm.add_constant(Z)

	results = sm.OLS(endog=yy, exog=Z).fit()

	beta_hat = results.params[2]
	gamma_hat = results.params[3]

	next_breakpoint_1 = current_breakpoint_1 - gamma_hat/beta_hat
	
	return next_breakpoint_1, results.params


def breakpoint_iterate():


	zz, yy = generate_breakpoint_data()
	current_breakpoint_1 = 1
	for i in range(6):
		print(i, current_breakpoint_1)
		current_breakpoint_1, params = breakpoint_fit(zz, yy, current_breakpoint_1)


	intercept = params[0]
	alpha_hat = params[1]
	beta_hat = params[2]
	breakpoint_hat = current_breakpoint_1

	yy_hats = intercept + alpha_hat*zz + beta_hat * np.maximum(zz - breakpoint_hat, 0)

	plt.plot(zz, yy_hats, linewidth=2, color="purple", linestyle="dashed")
	plt.scatter(zz, yy, s=4, color="green")

	plt.show()


breakpoint_iterate()