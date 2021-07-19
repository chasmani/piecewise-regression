
import numpy as np
import scipy
import math

def get_test_statistic(xx_davies, yy_davies, theta):
	"""
	Compute a test statistic for the Davies test for the p-value of existence of a breakpoint
	Based on Davies(1987)
		"Hypothesis Testing when a nuisance parameter is present only under the alternative"
	All the values worked out here are as named and described in that paper
	"""
	n = len(xx_davies)
	s_0 = 0
	s_1 = 0
	s_2 = 0
	s_3 = 0
	s_4 = 0

	for x in xx_davies:
		s_0 += x**2

		if x > theta:
			s_1 += x*(x-theta)
			s_3 += x - theta
	
		elif x < theta:
			s_2 += x*(x-theta)
			s_4 += x - theta

	a_hat = np.sum(yy_davies)/n
	b_hat = np.sum(xx_davies*yy_davies)/s_0
	V = s_1*s_2/s_0 + s_3*s_4/n
	S = 0
	for i in range(n):
		if xx_davies[i] > theta:
			S += (yy_davies[i] - a_hat - b_hat*xx_davies[i])*(xx_davies[i] - theta)
	S = S/(np.sqrt(np.abs(V)))
	return S


def davies_test(xx, yy, k=10, alternative="two_sided"):
	"""
	Significance test for the existence of a breakpoint
	Null hypothesis is that there is no breakpoint, or that the change in gradient is zero
	Alternative hypothesis is that there is a breakpoint, with a non-zero change in gradient
	The change is gradietn is a function of the breakpoint position
	The breakpoint posiition is a nuisannce parameter that only exists in the alternative hypothesis
	Based on Davies (1987), "Hypothesis Testing when a nuisance parameter is present only under the alternative"
	"""
	# Centre the x values - makes no difference to existence of a breakpoint
	# The Davies test has this as an assumption
	xx_davies = xx - np.mean(xx)
	yy_davies = yy
	
	# As in Muggeo's R package "segmented", cut from second to second to last data point
	# Need more data in the xx than in [L,U] for the test to work 
	L = xx_davies[1]
	U = xx_davies[-2]

	# More thetas is better
	thetas = np.linspace(L, U, k)
	# For each value of theta, compute a test statistic
	test_stats = []
	for theta in thetas:
		test_stat = get_test_statistic(xx_davies, yy_davies, theta)
		test_stats.append(test_stat)
	print(test_stats)
	if alternative == "two_sided":
		# Two sided test, M as defined by Davies
		M = np.max(np.abs(test_stats))
	elif alternative == "less":
		M = np.min(test_stats)
	elif alternative == "greater":
		M = np.max(test_stats)

	print("M is ", M)

	# Use formulas from Davies
	V = 0 
	for i in range(len(thetas)-1):
		V += np.abs(test_stats[i+1]-test_stats[i])

	print("V is ", V)
	p = scipy.stats.norm.cdf(-M) + V * np.exp(-.5*M**2) * 1/(np.sqrt(8*math.pi))
	
	print("p is ", p)

	if alternative == "two_sided":
		return p*2
	else:
		return p