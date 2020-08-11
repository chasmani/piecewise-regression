
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt

from probability_distributions import get_probabilities_power_law_finite_event_set


def convert_observations_into_ranked_empirical_counts(samples):

	counts = Counter(samples)
	n = [v for k,v in counts.most_common()]
	return n

#########################
# Fininte event set

def generate_samples_from_finite_power_law(exponent, W, N):
	"""
	generate random samples
	"""
	probs = get_probabilities_power_law_finite_event_set(exponent, W)
	return np.random.choice(W, size=N, p=probs)


def get_ranked_empirical_counts_from_finite_power_law(exponent, W, N):

	samples = generate_samples_from_finite_power_law(exponent, W, N)
	n = convert_observations_into_ranked_empirical_counts(samples)
	return n

######################
# Infinite event set

def generate_samples_from_infinite_power_law(exponent, N):

	xs_np = np.random.zipf(a=exponent, size = N)
	return xs_np


def get_ranked_empirical_counts_from_infinite_power_law(exponent, N):

	xs = generate_samples_from_infinite_power_law(exponent, N)
	n = convert_observations_into_ranked_empirical_counts(xs)
	return n





if __name__=="__main__":
	generate_samples_from_infinite_power_law(1.04, 200000)