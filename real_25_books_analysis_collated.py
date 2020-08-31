

import zipfanalysis

import os
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt

import re
import unidecode
from collections import Counter
from general_utilities import append_to_csv


def plot_books():

	directory = "gutenberg_manual"

	# 1. Get all files
	files = [f for f in listdir(directory)]
	print(files)
	plot_count = 1
	plt.figure(figsize=(16,12))
	for file in files:
		print(file)
		file_dir = "gutenberg_manual/" + file
		bookname = file.replace("_", " ").title()
		n = zipfanalysis.preprocessing.preprocessing.get_rank_frequency_from_text(file_dir)
		print(n)

		plt.subplot(5, 5, plot_count)

		book_name = file.replace(".txt","").replace("_"," ").title()

		plt.scatter(range(1, len(n)+1), n, s=1, label=book_name)
		plt.xscale("log")
		plt.yscale("log")
		#plt.xticks([])
		#plt.yticks([])

		plt.legend()

		if plot_count == 11:
			plt.ylabel("Frequency")

		if plot_count == 23:
			plt.xlabel("rank")

		plot_count += 1

	#plt.tight_layout()

	plt.savefig("images/top_25_books_rank_frequency.png")

	plt.show()

def plot_books_collated():

	directory = "gutenberg_manual"

	# 1. Get all files
	files = [f for f in listdir(directory)]
	print(files)
	plot_count = 1
	all_word_counts = Counter()
	for file in files:
		print(file)
		file_dir = "gutenberg_manual/" + file
		bookname = file.replace("_", " ").title()

		this_counts = zipfanalysis.preprocessing.preprocessing.get_word_counts(file_dir)
		all_word_counts += this_counts

	n = [v for k,v in all_word_counts.most_common()]
	plt.scatter(range(1, len(n)+1), n, s=2)
	plt.xscale("log")
	plt.yscale("log")

	plt.xlabel("rank")
	plt.ylabel("frequency")

	plt.savefig("images/top_25_books_collated_rank_frequency.png")
	plt.show()

def calculate_parameters_books():

	directory = "gutenberg_manual"
	results_file = "results/gutenberg_abc_analysis.csv"

	# 1. Get all files
	files = [f for f in listdir(directory)]
	print(files)
	plot_count = 1
	all_word_counts = Counter()
	for file in files:
		file_dir = "gutenberg_manual/" + file
		abc = zipfanalysis.abc(file_dir)
		n = zipfanalysis.preprocessing.preprocessing.get_rank_frequency_from_text(file_dir)
		V = len(n)
		N = sum(n)
		print(file, abc, V, N)
		csv_list = [file, abc, V, N]
		append_to_csv(csv_list, results_file)








	# 2. Remove project gutenberg bits

if __name__=="__main__":
	plot_books()