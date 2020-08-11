

import numpy as np
from scipy.special import zeta


def get_probabilities_power_law_finite_event_set(exponent, W):
	"""
	Generate a discrete power law probability distribution with exponent and W events
	"""
	probs = []

	for rank in range(1, W+1):
		prob = rank**(-1*exponent)
		probs.append(prob)
	normed_probs = np.array(probs)/sum(probs)
	return normed_probs


def get_probabilities_zeta_power_law(exponent, max_x):
	"""
	Generate the probabilites from a zeta distrbuted power law with inifinite event set, up to a maximum value
	"""
	probs = []
	z = zeta(exponent)
	for rank in range(1, max_x+1):
		prob = rank**(-1*exponent)/z
		probs.append(prob)
	return probs


if __name__=="__main__":
	p_1 = get_probabilities_power_law_finite_event_set(0.4, 6)
	print(p_1)	

	p_2 = get_probabilities_power_law_finite_event_set(-0.4, 6)
	print(p_2)