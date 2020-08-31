
import os
import re

import matplotlib.pyplot as plt
import pandas as pd
import powerlaw
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

def calculate_all_languages():

	for language_code in [
		#"chi-sim", 
		#"eng-gb",
		#"eng-fiction",
		"eng-all", 
		"eng-us-all", 
		"fre-all", 
		"ger-all", 
		"heb-all", 
		"ita-all", 
		"rus-all", 
		"spa-all"]:
		calculate_zipf(language_code)

def calculate_frequency_count_plaw_xmin_and_alpha_with_plot(language_code="eng-gb", year=1887):

	df = pd.read_csv("ngrams_by_year/{}-{}.csv".format(language_code, year), delimiter=";", names=["ngram", "year", "match_count", "volume_count", "language"])

	google_general_forms = [
		"_NOUN_",
		"_VERB_",
		"_ADP_",
		"_DET_",
		"_ADJ_",
		"_ADV_",
		"_PRON_",
		"_NUM_",
		"_CONJ_",
		"_PRT_"
	]


	print(df)

	df = df[~df['ngram'].astype(str).str.contains('\_')]

	# Only look at words that appeared 10 times or more - to make it possible to compute
	

	# Only look at top 30000 words
	df_top = df.nlargest(50000, "match_count").reset_index(drop=True).reset_index()

	print("Now: ", df_top)

	df = df.sort_values(by="match_count", ascending=False).reset_index(drop=True).reset_index()

	print(df)


	zs = df_top["match_count"].to_numpy()

	# The way powerlaw library works, the zs are the data format to approximate the f(z) power law
	lib_fit = powerlaw.Fit(zs, discrete=True, estimate_discrete=True)
	alpha = lib_fit.alpha
	fmin = int(lib_fit.xmin)
	# Maximum ranked word in the powerlaw regime
	#r_max = sum(f[xmin-1:])
	
	print(alpha, fmin)



	df = df.sort_values(by="match_count", ascending=False).reset_index(drop=True).reset_index()

	df_c = df.groupby(by="match_count").count().reset_index()
	z_axis = df_c["match_count"].to_numpy()
	f = df_c["ngram"].to_numpy()
	
	# Get df of frequency counts within fmin range
	# This df includes freqnecuies and number of tokens with that frequency
	# TOtal tokens in the range is the sum of the counts within the range
	df_region = df_c[df_c["match_count"] >= fmin]
	rmax = df_region["ngram"].sum()
	print(rmax)

	# For the plot
	df_region["f_hats"] = df_region['match_count'].apply(lambda x: 10**11*x**(-1*alpha))
	
	print(df_region)

	plt.scatter(x=z_axis, y=f, s=5)

	plt.xscale("log")
	plt.ylim(0.8, 1000000)

	plt.yscale("log")

	plt.xlabel("number of occurences of token, z")
	plt.ylabel("$\phi(z)$ - Number of tokens with this count")
	plt.title("Counts Distribution of Google 1-grams {}".format(year))

	plt.plot(df_region["match_count"], df_region["f_hats"], color="red")

	plt.savefig("images/year_words_only_count_dist_plot_{}_{}.png".format(language_code, year))
	
	plt.show()
	

def calculate_frequency_count_plaw_xmin_and_alpha(language_code="eng-fiction", filename=""):

	year = re.findall(r'\d+', filename)[0]
	print(filename, year)

	df = pd.read_csv("ngrams_by_year/{}-{}.csv".format(language_code, year), delimiter=";", names=["ngram", "year", "match_count", "volume_count", "language"])

	google_general_forms = [
		"_NOUN_",
		"_VERB_",
		"_ADP_",
		"_DET_",
		"_ADJ_",
		"_ADV_",
		"_PRON_",
		"_NUM_",
		"_CONJ_",
		"_PRT_"
	]


	df = df[~df['ngram'].astype(str).str.contains('\_')]
	# Only look at top 30000 words
	df_top = df.nlargest(30000, "match_count").reset_index(drop=True).reset_index()

	zs = df_top["match_count"].to_numpy()

	# The way powerlaw library works, the zs are the data format to approximate the f(z) power law
	lib_fit = powerlaw.Fit(zs, discrete=True, estimate_discrete=True)
	alpha = lib_fit.alpha
	fmin = int(lib_fit.xmin)

	df_sorted = df_top.sort_values(by="match_count", ascending=False).reset_index(drop=True).reset_index()

	df_c = df_sorted.groupby(by="match_count").count().reset_index()
	z_axis = df_c["match_count"].to_numpy()
	f = df_c["ngram"].to_numpy()
	
	# Get df of frequency counts within fmin range
	# This df includes freqnecuies and number of tokens with that frequency
	# Total tokens in the range is the sum of the counts within the range
	df_region = df_c[df_c["match_count"] >= fmin]
	rmax = df_region["ngram"].sum()
	print(alpha, fmin, rmax)

	csv_row = [language_code, year, "Better Try 30000 top words. Powerlaw beta fmin r_max on f(z)", alpha, fmin, rmax]
	append_to_csv(csv_row, "results/powerlaw_zipf_results_by_year.csv")


def calculate_powerlaw_lib(language_code="eng-gb"):
	
	ngrams_by_year_files = os.listdir("ngrams_by_year")
	this_language_files = [f for f in ngrams_by_year_files if language_code in f]

	print(this_language_files)

	for yearfile in sorted(this_language_files):
		calculate_frequency_count_plaw_xmin_and_alpha(language_code, "ngrams_by_year/"+yearfile)

def powerlib_all():

	for language_code in [
		"chi-sim", 
		"eng-gb",
		"eng-fiction",
		"eng-all", 
		"eng-us-all", 
		"fre-all", 
		"ger-all", 
		"heb-all", 
		"ita-all", 
		"rus-all", 
		"spa-all"]:
		calculate_powerlaw_lib(language_code)


if __name__=="__main__":
	powerlib_all()
