
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import math


def test_it():

	alpha = -4
	beta_1 = -2
	intercept = 100
	breakpoint_1 = 7

	n_points = 200

	xx = np.linspace(0, 20, n_points)
	yy = intercept + alpha*xx + beta_1 * np.maximum(xx - breakpoint_1, 0) + np.random.normal(size=n_points)

def test_null_basic():

	

	intercept = 5
	alpha = 1
	beta_1 = 0
	breakpoint_1 = 2
	n_points = 200

	xx = np.linspace(-9.5, 9.5, n_points)
	yy = intercept + alpha*xx + beta_1 * np.maximum(xx - breakpoint_1, 0) + np.random.normal(size=n_points)

	L = -8
	U = 8 


	thetas = np.arange(L, U, 0.2)

	test_stats = []

	for theta in thetas:
		test_stat = get_test_statistic(xx, yy, theta)
		test_stats.append(test_stat)

		# Two sided
	M = np.max(np.abs(test_stats))

		# For one sided test - also don't multiply p by 2
	#M = np.max(test_stats)

	V = 0 
	for i in range(len(thetas)-1):
		V += np.abs(test_stats[i+1]-test_stats[i])

	p = norm.cdf(-M) + V * np.exp(-.5*M**2) * 1/(np.sqrt(8*math.pi))
	p = p*2

	return p



def get_test_statistic(xx, yy, theta):

	n = len(xx)

	s_0 = 0
	s_1 = 0
	s_2 = 0
	s_3 = 0
	s_4 = 0

	for x in xx:
		s_0 += x**2

		if x > theta:
			s_1 += x*(x-theta)
			s_3 += x - theta

		elif x < theta:
			s_2 += x*(x-theta)
			s_4 += x - theta


	a_hat = np.sum(yy)/n
	b_hat = np.sum(xx*yy)/s_0

	#print(a_hat, b_hat)

	V = s_1*s_2/s_0 + s_3*s_4/n

	S = 0
	for i in range(n):
		if xx[i] > theta:
			S += (yy[i] - a_hat - b_hat*xx[i])*(xx[i] - theta)

	S = S/(np.sqrt(V))

	return S

def check_p_values():

	p_count_5 = 0 
	p_count_2 = 0
	p_count_1 = 0

	for seed in range(1000):
		np.random.seed(seed)

		p = test_null_basic()

		if p<0.05:
			p_count_5 +=1
		if p<0.02:
				p_count_2 += 1
		if p<0.01:
			p_count_1 +=1

	print("P values for 5, 2, and 1 percent are ")
	print(p_count_5/10, p_count_2/10, p_count_1/10)

	#plt.scatter(thetas, test_stats)
	#plt.show()




check_p_values()