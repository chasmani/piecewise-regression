
import os
import re

import matplotlib.pyplot as plt
import pandas as pd
import zipfanalysis

from general_utilities import append_to_csv


def calculate_zipf(language_code="eng-gb"):
	
	ngrams_by_year_files = os.listdir("ngrams_by_year")
	this_language_files = [f for f in ngrams_by_year_files if language_code in f]

	print(this_language_files)

	for yearfile in sorted(this_language_files):
		calculate_zipf_from_file(language_code, "ngrams_by_year/"+yearfile)

def plot_zipf_from_file(language, filename):

		year = re.findall(r'\d+', filename)[0]
		print(filename, year)

		df = pd.read_csv(filename, delimiter=";", names=["ngram", "year", "match_count", "volume_count", "language"])

		# Remove duplicates
		df = df.drop_duplicates()

		# Remove any strings you don't want
		df = df[df['ngram'] != ","]


		df = df.sort_values(by="match_count", ascending=False).reset_index(drop=True).reset_index()

		df["rank"] = df["index"] + 1

		print(df.head())

		df.plot.scatter(x="rank", y="match_count", s=5)

		plt.xscale("log")
		plt.yscale("log")

		plt.xlabel("Rank")
		plt.ylabel("Frequency")
		plt.title("Rank-Frequency {} {}".format(language, year))
		plt.savefig("images/rank_frequency_{}_{}.png".format(language, year))


def calculate_zipf_from_file(language, filename):

		year = re.findall(r'\d+', filename)[0]
		print(filename, year)

		df = pd.read_csv(filename, delimiter=";", names=["ngram", "year", "match_count", "volume_count", "language"])

		# Remove duplicates
		df = df.drop_duplicates()

		# Remove any strings you don't want
		df = df[df['ngram'] != ","]

		df = df.nlargest(5000, "match_count")
		
		# Get n vector of largest 5000 counts
		ns = df["match_count"].astype(int).to_numpy()

		alpha_clauset = zipfanalysis.estimators.clauset.clauset_estimator(ns)
		csv_row = [language, year, "Clauset top 5000 types", alpha_clauset]
		append_to_csv(csv_row, "results/zipf_results_by_year.csv")

		alpha_pdf = zipfanalysis.estimators.ols_regression_pdf.ols_regression_pdf_estimator(ns)
		csv_row = [language, year, "PDF Regression top 5000 types", alpha_pdf]
		append_to_csv(csv_row, "results/zipf_results_by_year.csv")


if __name__=="__main__":
	calculate_zipf("eng-fiction")
	calculate_zipf("eng-gb")