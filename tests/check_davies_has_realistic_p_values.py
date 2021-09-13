
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import math

import os, sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))

from piecewise_regression.davies import davies_test

def check_p_values():

	p_count_20 = 0
	p_count_5 = 0 
	p_count_2 = 0
	p_count_1 = 0

	p_values = []

	sample_size = 1000

	for seed in range(sample_size):
		np.random.seed(seed)

		xx_bp, yy_bp = generate_data()
		theta = 0

		p = davies_test(xx_bp, yy_bp)

		if p<0.2:
			p_count_20 += 1
		if p<0.05:
			p_count_5 +=1
		if p<0.02:
				p_count_2 += 1
		if p<0.01:
			p_count_1 +=1

		p_values.append(p)

	#plt.hist(p_values)
	#plt.show()


	print("P values from empirical testing, along with excpetd p-values")
	print("p-values \t\t\t0.2  \t0.05\t0.02\t0.01")
	print("fraction less than \t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}".format(p_count_20/sample_size, p_count_5/sample_size,
		p_count_2/sample_size,p_count_1/sample_size))

	print("The empirical fractions should be slghtly less than the expected p values")


def check_p_values_less():

	p_count_20 = 0
	p_count_5 = 0 
	p_count_2 = 0
	p_count_1 = 0

	p_values = []

	sample_size = 1000

	for seed in range(sample_size):
		np.random.seed(seed)

		xx_bp, yy_bp = generate_data()
		theta = 0

		p = davies_test(xx_bp, yy_bp, alternative="less")

		if p<0.2:
			p_count_20 += 1
		if p<0.05:
			p_count_5 +=1
		if p<0.02:
				p_count_2 += 1
		if p<0.01:
			p_count_1 +=1

		p_values.append(p)

	#plt.hist(p_values)
	#plt.show()


	print("P values from empirical testing, along with excpetd p-values")
	print("p-values \t\t\t0.2  \t0.05\t0.02\t0.01")
	print("fraction less than \t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}".format(p_count_20/sample_size, p_count_5/sample_size,
		p_count_2/sample_size,p_count_1/sample_size))

	print("The empirical fractions should be slightly less than the expected p values")


def generate_data():

	intercept = 5
	alpha = 1
	beta_1 = 0
	breakpoint_1 = 2
	n_points = 20

	xx_bp_1 = list(np.linspace(-9.5, 9.5, n_points))
	xx_bp_2 = list(np.linspace(0,9.5,n_points))
	xx_bp = np.array(xx_bp_1 + xx_bp_2)
	#xx_bp = np.array(xx_bp_1)

	yy_bp = intercept + alpha*xx_bp + beta_1 * np.maximum(xx_bp - breakpoint_1, 0) + np.random.normal(size=len(xx_bp))

	return xx_bp, yy_bp


if __name__=="__main__":
	check_p_values()
	check_p_values_less()