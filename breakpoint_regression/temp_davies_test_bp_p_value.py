
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

	xx_bp = np.linspace(0, 20, n_points)
	yy_bp = intercept + alpha*xx_bp + beta_1 * np.maximum(xx_bp - breakpoint_1, 0) + np.random.normal(size=n_points)

def test_null_basic_centered():

	intercept = 5
	alpha = 1
	beta_1 = 0
	breakpoint_1 = 2
	n_points = 200

	xx_bp = np.linspace(-9.5, 9.5, n_points)
	yy_bp = intercept + alpha*xx_bp + beta_1 * np.maximum(xx_bp - breakpoint_1, 0) + np.random.normal(size=n_points)

	L = -8
	U = 8 


	thetas = np.arange(L, U, 0.2)

	test_stats = []

	for theta in thetas:
		test_stat = get_test_statistic(xx_bp, yy_bp, theta)
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



def get_test_statistic(xx_bp, yy_bp, theta):

	n = len(xx_bp)

	s_0 = 0
	s_1 = 0
	s_2 = 0
	s_3 = 0
	s_4 = 0

	for x in xx_bp:
		s_0 += x**2

		if x > theta:
			s_1 += x*(x-theta)
			s_3 += x - theta

		elif x < theta:
			s_2 += x*(x-theta)
			s_4 += x - theta


	a_hat = np.sum(yy_bp)/n
	b_hat = np.sum(xx_bp*yy_bp)/s_0

	#print(a_hat, b_hat)

	#print(s_0, s_1, s_2, s_3, s_4)

	V = s_1*s_2/s_0 + s_3*s_4/n

	#print("V is ", V)

	S = 0
	for i in range(n):
		if xx_bp[i] > theta:
			S += (yy_bp[i] - a_hat - b_hat*xx_bp[i])*(xx_bp[i] - theta)

	S = S/(np.sqrt(np.abs(V)))

	return S

def check_p_values():

	p_count_5 = 0 
	p_count_2 = 0
	p_count_1 = 0

	p_values = []

	for seed in range(100):
		np.random.seed(seed)

		xx_bp, yy_bp = generate_data()
		L = np.percentile(xx_bp, 10)
		U = np.percentile(xx_bp, 90)
		theta = 0

		p = davies_test_2_sided(xx_bp, yy_bp)
		print(p)

		if p<0.05:
			p_count_5 +=1
		if p<0.02:
				p_count_2 += 1
		if p<0.01:
			p_count_1 +=1

		p_values.append(p)

	#plt.hist(p_values)
	#plt.show()


	print("P values for 5, 2, and 1 percent are ")
	print(p_count_5, p_count_2, p_count_1)

	#plt.scatter(thetas, test_stats)
	#plt.show()

def generate_data():

	intercept = 5
	alpha = 1
	beta_1 = 0
	breakpoint_1 = 2
	n_points = 7

	xx_bp_1 = list(np.linspace(-9.5, 9.5, n_points))
	xx_bp_2 = list(np.linspace(0,9.5,n_points))
	xx_bp = np.array(xx_bp_1 + xx_bp_2)
	#xx_bp = np.array(xx_bp_1)

	yy_bp = intercept + alpha*xx_bp + beta_1 * np.maximum(xx_bp - breakpoint_1, 0) + np.random.normal(size=len(xx_bp))

	return xx_bp, yy_bp


def davies_test_2_sided(xx_bp, yy_bp):
	"""
	Significane test for the existence of a breakpoint
	Null hypothesis is that there is no breakpoint, or that the change in gradient is zero
	Alternative hypothesis is that there is a breakpoint, with a non-zero change in gradient
	The change is gradietn is a function of the breakpoint position
	The breakpoint posiition is a nuisannce parameter that only exists in the alternative hypothesis
	Based on Davies (1987), "Hypothesis Testing when a nuisance parameter is present only under the alternative"
	"""

	# Centre the x values - makes no difference to existence of a breakpoint
	# The Davies test has this as an assumption
	xx_bp = xx_bp - np.mean(xx_bp)
	
	# Cut into the data a little bit to avoid errors caused by the constants being zero
	# Go a bit wider (x1.11) on selecting the region to counteract this.
	# You need a data region slightly bigger than the [L,U] region for this to work 
	L = np.percentile(xx_bp, 10)
	U = np.percentile(xx_bp, 90)

	# More thetas is better
	thetas = np.linspace(L, U, 20)

	# For each value of theta, compute a test statistic
	test_stats = []
	for theta in thetas:
		test_stat = get_test_statistic(xx_bp, yy_bp, theta)
		test_stats.append(test_stat)

	# Two sided test, M as defined by Davies
	M = np.max(np.abs(test_stats))

	# For one sided test - also don't multiply p by 2
	#M = np.max(test_stats)

	# Use formulas from Davies and Muggeo
	V = 0 
	for i in range(len(thetas)-1):
		V += np.abs(test_stats[i+1]-test_stats[i])

	p = norm.cdf(-M) + V * np.exp(-.5*M**2) * 1/(np.sqrt(8*math.pi))
	# Two sided test, beta can be positive or negative
	p = p*2

	return p







check_p_values()