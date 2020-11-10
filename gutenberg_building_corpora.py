
import pandas as pd
import zipfile
import matplotlib.pyplot as plt
import csv
from io import TextIOWrapper
from collections import Counter
import numpy as np
import seaborn as sns

pd.set_option('display.max_columns', None)


def build_corpora():

	np.random.seed(1)

	n_books = 3

	input_file = "results/gutenberg_analysis_no_errors.csv"

	df = pd.read_csv(input_file, sep=";", names=["gut_id", "title", "author", "author_birth_year", 
			"author_birth_death", "language", "downloads", "subjects", "filetype", "alpha", "total_words", "unique_words"])

	df["alpha"] = pd.to_numeric(df["alpha"], errors="coerce")
	df["unique_words"] = pd.to_numeric(df["unique_words"], errors="coerce")
	df["total_words"] = pd.to_numeric(df["total_words"], errors="coerce")

	

	
	df = df[df["language"] == "['en']"]
	"""
	df = df[df["total_words"] > 1000]
	"""

	sample_df = df.sample(n_books)
	print(sample_df)

	sample_df["alpha"].hist()
	#plt.show()



	# Prediction 1 
	# Alpha will be lower than the mean of alphas

	# Prediction 2
	# Breakpoint will be approximately N_avg

	counters = []
	corpus_counter = Counter()

	df = pd.DataFrame(columns = ["ngram", "count"])

	for uid, gut_id in sample_df["gut_id"].iteritems():
		
		this_counter = Counter()

		zf = zipfile.ZipFile('gutenberg_standard/SPGC-counts-2018-07-18.zip') 
		this_df = pd.read_csv(zf.open('SPGC-counts-2018-07-18/{}_counts.txt'.format(gut_id)), delimiter="\t", names=["ngram", "count"])
		
		df = pd.merge(df, this_df, on="ngram", how="outer", suffixes=[None, gut_id])


	df["count"] = df.sum(axis=1)





	print(df)


def open_gutenberg_file_csv_reader(gut_id):

	with zipfile.ZipFile('gutenberg_standard/SPGC-counts-2018-07-18.zip') as zf:
		with zf.open('SPGC-counts-2018-07-18/{}_counts.txt'.format(gut_id)) as infile:
			reader = csv.reader(TextIOWrapper(infile, "utf-8"), delimiter="\t")
			for row in reader:
				print(row)


def entropy_sums():

	np.random.seed(1)

	n_books = 3

	input_file = "results/gutenberg_analysis_no_errors.csv"

	df = pd.read_csv(input_file, sep=";", names=["gut_id", "title", "author", "author_birth_year", 
			"author_birth_death", "language", "downloads", "subjects", "filetype", "alpha", "total_words", "unique_words"])

	df["alpha"] = pd.to_numeric(df["alpha"], errors="coerce")
	df["unique_words"] = pd.to_numeric(df["unique_words"], errors="coerce")
	df["total_words"] = pd.to_numeric(df["total_words"], errors="coerce")

	

	
	df = df[df["language"] == "['en']"]
	"""
	df = df[df["total_words"] > 1000]
	"""

	sample_df = df.sample(n_books)
	print(sample_df)

	#sample_df["alpha"].hist()
	#plt.show()



	# Prediction 1 
	# Alpha will be lower than the mean of alphas

	# Prediction 2
	# Breakpoint will be approximately N_avg

	counters = []
	corpus_counter = Counter()

	df = pd.DataFrame(columns = ["ngram", "count"])

	for uid, gut_id in sample_df["gut_id"].iteritems():
		
		this_counter = Counter()

		zf = zipfile.ZipFile('gutenberg_standard/SPGC-counts-2018-07-18.zip') 
		this_df = pd.read_csv(zf.open('SPGC-counts-2018-07-18/{}_counts.txt'.format(gut_id)), delimiter="\t", names=["ngram", "count"])

		count_sum = this_df["count"].sum()
		#this_df["normalised_freq"] = this_df["count"]/count_sum
		#this_df["entropy"] = this_df["normalised_freq"].apply(lambda x: x*np.log(x))

		print(this_df)

		
		df = pd.merge(df, this_df, on="ngram", how="outer", suffixes=[None, gut_id])


	df["count"] = df.sum(axis=1)
	df = df.sort_values(by="count", ascending=False).reset_index(drop=True).reset_index()
	df["rank"] = df["index"] + 1
	gut_ids = [""] + list(sample_df["gut_id"].to_numpy())

	print(gut_ids)

	for gut_id in gut_ids:
		print("Gut id is ", gut_id)
		this_book_sum = df["count{}".format(gut_id)].sum()
		df["normalised_freq{}".format(gut_id)] = df["count{}".format(gut_id)]/this_book_sum
		df["entropy{}".format(gut_id)] = df["normalised_freq{}".format(gut_id)].apply(lambda x: -x*np.log(x))
		df["total_entropy{}".format(gut_id)] = df["entropy{}".format(gut_id)].cumsum()
		df["cum_words{}".format(gut_id)] = df["count{}".format(gut_id)].cumsum()


	print(df)

	for gut_id in gut_ids:

		if gut_id == "":
			label="Corpora"
		else:
			label=gut_id

		sns.lineplot(x="rank", y="total_entropy{}".format(gut_id), data=df, label=label)

	plt.xlabel("Corpus Rank")
	plt.ylabel("Total Entropy")

	plt.xscale("log")


	plt.legend()

	plt.show()


def entropy_sums(seed, n_books):

	np.random.seed(seed)



	input_file = "results/gutenberg_analysis_no_errors.csv"

	df = pd.read_csv(input_file, sep=";", names=["gut_id", "title", "author", "author_birth_year", 
			"author_birth_death", "language", "downloads", "subjects", "filetype", "alpha", "total_words", "unique_words"])

	df["alpha"] = pd.to_numeric(df["alpha"], errors="coerce")
	df["unique_words"] = pd.to_numeric(df["unique_words"], errors="coerce")
	df["total_words"] = pd.to_numeric(df["total_words"], errors="coerce")

	

	
	df = df[df["language"] == "['en']"]
	"""
	df = df[df["total_words"] > 1000]
	"""

	sample_df = df.sample(n_books)


	#sample_df["alpha"].hist()
	#plt.show()



	# Prediction 1 
	# Alpha will be lower than the mean of alphas

	# Prediction 2
	# Breakpoint will be approximately N_avg

	counters = []
	corpus_counter = Counter()

	df = pd.DataFrame(columns = ["ngram", "count"])

	for uid, gut_id in sample_df["gut_id"].iteritems():
		
		this_counter = Counter()

		zf = zipfile.ZipFile('gutenberg_standard/SPGC-counts-2018-07-18.zip') 
		this_df = pd.read_csv(zf.open('SPGC-counts-2018-07-18/{}_counts.txt'.format(gut_id)), delimiter="\t", names=["ngram", "count"])

		count_sum = this_df["count"].sum()
		#this_df["normalised_freq"] = this_df["count"]/count_sum
		#this_df["entropy"] = this_df["normalised_freq"].apply(lambda x: x*np.log(x))


		
		df = pd.merge(df, this_df, on="ngram", how="outer", suffixes=[None, gut_id])


	df["count"] = df.sum(axis=1)
	df = df.sort_values(by="count", ascending=False).reset_index(drop=True).reset_index()
	df["rank"] = df["index"] + 1
	gut_ids = [""] + list(sample_df["gut_id"].to_numpy())


	for gut_id in gut_ids:

		this_book_sum = df["count{}".format(gut_id)].sum()
		df["normalised_freq{}".format(gut_id)] = df["count{}".format(gut_id)]/this_book_sum
		df["entropy{}".format(gut_id)] = df["normalised_freq{}".format(gut_id)].apply(lambda x: -x*np.log(x))
		df["total_entropy{}".format(gut_id)] = df["entropy{}".format(gut_id)].cumsum()
		df["cum_words{}".format(gut_id)] = df["count{}".format(gut_id)].cumsum()

	max_rank = df["rank"].max()


	all_Hs = []
	premix_Hs = []


	counts = df["count"].to_numpy()

	ranks = range(1, max_rank+1, 1000)

	gut_ids_books_only = list(sample_df["gut_id"].to_numpy())

	for rank in ranks:


		all_counts = counts[:rank]

		all_H = get_entropy_of_counts(all_counts)
		all_Hs.append(all_H)

		weighted_H = 0

		for gut_id in gut_ids_books_only:
			gut_counts = df["count{}".format(gut_id)].to_numpy()[:rank]

			gut_counts = np.nan_to_num(gut_counts, nan=0)

			gut_H = get_entropy_of_counts(gut_counts)

			weighted_H += sum(gut_counts)*gut_H
		premix_H = weighted_H / sum(all_counts)


		premix_Hs.append(premix_H)


	plt.plot(ranks, all_Hs, label="Post Mix H")
	plt.plot(ranks, premix_Hs, label="Pre Mix H")
	plt.xscale("log")
	plt.xlabel("Corpus Rank")
	plt.ylabel("Entropy Up To Rank")
	plt.legend()
	plt.title("Entropy of Corpus Increases During Composition\n{} books randomly selected with seed {}".format(n_books, seed))
	plt.savefig("images/entropy_pre_and_post_mix_seed_{}_n_{}.png".format(seed, n_books))

	plt.close()


def get_entropy_of_counts(ns):

	N = sum(ns)

	H = 0
	for n in ns:
		if not np.isnan(n) and n>0:
			H += -n/N * np.log(n/N)
	return H

def plot_loads_pre_post_mix():

	for n_books in [1,2,5,10,20,50,100,200]:
		for seed in [1,2,3,4,5]:
			print(n_books, seed)
			entropy_sums(seed, n_books)



plot_loads_pre_post_mix()