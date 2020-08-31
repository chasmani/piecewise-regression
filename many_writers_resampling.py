
from collections import Counter
import csv

import numpy as np
import matplotlib.pyplot as plt

from data_generators import generate_samples_from_infinite_power_law, convert_observations_into_ranked_empirical_counts
from probability_distributions import get_probabilities_power_law_finite_event_set


def one_writer():

	words = generate_samples_from_infinite_power_law(1.1, 100000)
	print(words)
	word_counter = Counter(words)


	ranked = convert_observations_into_ranked_empirical_counts(word_counter)

	print(ranked)

	plt.scatter(range(1, len(ranked)+1), ranked)


	plt.xscale("log")
	plt.yscale("log")


	plt.show()


def one_writer_words_and_probs(exponent, W):

	words = generate_samples_from_infinite_power_law(exponent, W)
	
	word_counter = Counter(words)

	words = [k for k,v in word_counter.most_common()]
	counts = [v for k,v in word_counter.most_common()]
	probs = np.array(counts)/sum(counts)




	print(words, counts)

	return words, probs


def resample_from_one_writer(exponent, W, N, number_of_books):

	words, probs = one_writer_words_and_probs(exponent, W)

	samples = []
	for i in range(number_of_books):

		sample = np.random.choice(words, p=probs, size=N)
		samples.append(sample)


	return samples


def plot_counter_rank_freq(counter):
	ranked = [v for k,v in counter.most_common()]
	plt.scatter(range(1, len(ranked)+1), ranked)
	plt.xscale("log")
	plt.yscale("log")

	plt.xlabel("Rank")
	plt.ylabel("Frequency")

	plt.show()



def lots_of_samples():

	full_count = Counter()

	all_writer_samples = []
	for k in range(1):
		writer_samples = resample_from_one_writer(1.1, 10000, 1000, 4)
		all_writer_samples += writer_samples

	for sample in all_writer_samples:
		full_count.update(sample)
	

	plot_counter_rank_freq(full_count)


def get_writer_words_and_power_law_probs(exponent = 1.1, W=10000):

	words = generate_samples_from_infinite_power_law(exponent, W)


	word_counter = Counter(words)
	words = [k for k,v in word_counter.most_common()]

	print(len(words))

	
	probs = get_probabilities_power_law_finite_event_set(exponent, len(words))

	return words, probs





def lots_of_samples_with_individual_power_laws():

	full_count = Counter()
	exponent = 1.2
	W = 100000
	N = 50000
	writer_count = 25

	for writer in range(writer_count):

		words, probs = get_writer_words_and_power_law_probs(exponent, W=W)

		for books in range(1):
			sample = np.random.choice(words, p=probs, size=N)
			full_count.update(sample)


	ranked = [v for k,v in full_count.most_common()]
	
	plt.scatter(range(1, len(ranked)+1), ranked, s=10)
	plt.xscale("log")
	plt.yscale("log")

	plt.xlabel("Rank")
	plt.ylabel("Frequency")

	plt.title("Zipf's Law {} Writers\ninitial W={}, exponent={}, book size={}".format(writer_count, W, exponent, N))

	plt.savefig("images/zipf_law_kinked_{}_authors_W_{}_lamb_{}_N_{}.png".format(writer_count, W, exponent, N))
	plt.show()



def frequency_count_representation():


	full_count = Counter()
	exponent = 1.2
	W = 100000
	N = 100000
	writer_count = 50

	for writer in range(writer_count):

		words, probs = get_writer_words_and_power_law_probs(exponent, W=W)

		for books in range(1):
			sample = np.random.choice(words, p=probs, size=N)
			full_count.update(sample)


	ranked = [v for k,v in full_count.most_common()]
	print(ranked)

	f_n = Counter(ranked)
	print(f_n)

	n = [k for k,v in f_n.most_common()]
	f_n = [v for k,v in f_n.most_common()]

	plt.scatter(n, f_n, s=10)

	plt.xscale("log")
	plt.yscale("log")


	plt.show()

def volume_vs_frequency_rank_sim():


	full_count = Counter()
	volume_count = Counter()
	exponent = 1.2
	W = 100000
	N = 100000
	writer_count = 50

	for writer in range(writer_count):

		words, probs = get_writer_words_and_power_law_probs(exponent, W=W)

		for books in range(1):
			sample = np.random.choice(words, p=probs, size=N)
			full_count.update(sample)
			unique_words = list(set(sample))
			volume_count.update(unique_words)

	volumes_of_ranks = []
	for k,v in full_count.most_common():
		volumes_of_ranks.append(volume_count[k])


	plt.scatter(range(1, len(volumes_of_ranks)+1), volumes_of_ranks, s=2)
	plt.xscale("log")
	plt.yscale("log")

	plt.xlabel("Frequency Rank")
	plt.ylabel("Volume Count")

	plt.title("Volume Count vs Frequency Rank {} Writers\ninitial W={}, exponent={}, book size={}".format(writer_count, W, exponent, N))

	plt.savefig("images/sim_volume_vs_freq_rankzipf_law_kinked_{}_authors_W_{}_lamb_{}_N_{}.png".format(writer_count, W, exponent, N))
	plt.show()


def lots_of_samples_based_on_top_25_books_subplots():

	W = 100000
	input_csv = "results/gutenberg_abc_analysis.csv"
	full_count = Counter()
	plot_count = 1
	plt.figure(figsize=(16,12))
	with open(input_csv) as csvfile:
		reader = csv.reader(csvfile, delimiter=';')
		for row in reader:
			alpha = float(row[1])
			N = int(row[3])

			
			words, probs = get_writer_words_and_power_law_probs(alpha, W=W)
			sample = np.random.choice(words, p=probs, size=N)
			word_counts = Counter(sample)
			
			n = [v for k,v in word_counts.most_common()]

			#full_count.update(sample)

			plt.subplot(5, 5, plot_count)

			book_name = row[0].replace(".txt","").replace("_"," ").title()

			plt.scatter(range(1, len(n)+1), n, s=1, label=book_name)
			plt.xscale("log")
			plt.yscale("log")

			plt.legend()

			if plot_count == 11:
				plt.ylabel("Frequency")

			if plot_count == 23:
				plt.xlabel("rank")

			plot_count += 1

	plt.savefig("images/top_25_books_rank_frequency_simulated.png")

	plt.show()

def lots_of_samples_based_on_top_25_books_collated():

	W = 200000
	input_csv = "results/gutenberg_abc_analysis.csv"
	full_count = Counter()
	plot_count = 1

	with open(input_csv) as csvfile:
		reader = csv.reader(csvfile, delimiter=';')
		for row in reader:
			alpha = float(row[1])
			N = int(row[3])

			
			words, probs = get_writer_words_and_power_law_probs(alpha, W=W)
			sample = np.random.choice(words, p=probs, size=N)
			


			full_count.update(sample)

	n = [v for k,v in full_count.most_common()]

			
	plt.scatter(range(1, len(n)+1), n, s=2)
	plt.xscale("log")
	plt.yscale("log")

	plt.legend()

	plt.ylabel("frequency")

	plt.xlabel("rank")

	plt.savefig("images/top_25_books_collated_simulatied_rank_frequency.png")

	plt.show()


if __name__=="__main__":
	lots_of_samples_based_on_top_25_books_subplots()