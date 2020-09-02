
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




def generate_double_breakpoint_data():

	alpha = -4
	beta_1 = -2
	beta_2 = -2
	intercept = 100
	breakpoint_1 = 7
	breakpoint_2 = 10

	n_points = 200

	xx = np.linspace(0, 20, n_points)

	yy = intercept + alpha*xx + beta_1 * np.maximum(xx - breakpoint_1, 0) + beta_2 * np.maximum(xx - breakpoint_2, 0)  + np.random.normal(size=n_points)

	return xx, yy

	
	#plt.scatter(xx, yy, s=2)
	#plt.show()



def breakpoint_fit(zz, yy, current_breakpoint_1=8, current_breakpoint_2=10.5):

	A = zz
	U1 = (zz - current_breakpoint_1) * np.heaviside(zz- current_breakpoint_1, 1) 
	V1 = np.heaviside(zz - current_breakpoint_1, 1)

	U2 = (zz - current_breakpoint_2) * np.heaviside(zz- current_breakpoint_2, 1)
	V2 = np.heaviside(zz - current_breakpoint_2, 1)

	
	Z = np.array([A, U1 , U2, V1, V2]).T

	Z = sm.add_constant(Z)

	results = sm.OLS(endog=yy, exog=Z).fit()

	
	beta_1_hat = results.params[2]
	beta_2_hat = results.params[3]
	gamma_1_hat = results.params[4]
	gamma_2_hat = results.params[5]
	
	next_breakpoint_1 = current_breakpoint_1 - gamma_1_hat/beta_1_hat
	
	next_breakpoint_2 = current_breakpoint_2 - gamma_2_hat/beta_2_hat
	
	return next_breakpoint_1, next_breakpoint_2, results.params
	

def double_breakpoint_iterate():


	zz, yy = generate_double_breakpoint_data()
	current_breakpoint_1 = 9
	current_breakpoint_2 = 9.1




	for i in range(6):
		current_breakpoint_1, current_breakpoint_2, params = breakpoint_fit(zz, yy, current_breakpoint_1, current_breakpoint_2)
		print(current_breakpoint_1, current_breakpoint_2)

	intercept = params[0]
	alpha_hat = params[1]
	beta_1_hat = params[2]
	beta_2_hat = params[3]
	breakpoint_1_hat = current_breakpoint_1
	breakpoint_2_hat = current_breakpoint_2

	yy_hats = intercept + alpha_hat*zz + beta_1_hat * np.maximum(zz - breakpoint_1_hat, 0) + beta_2_hat * np.maximum(zz - breakpoint_2_hat, 0)

	plt.plot(zz, yy_hats, linewidth=2, color="purple", linestyle="dashed")
	plt.scatter(zz, yy, s=4, color="green")

	plt.show()
	
	



double_breakpoint_iterate()