
import pandas as pd
import matplotlib.pyplot as plt

from probability_distributions import get_probabilities_zeta_power_law

def plot_all_years_exponents(min_year = 1800, language_code="eng-fiction"):

	results_filename = "results/powerlaw_zipf_results_by_year.csv"

	df = pd.read_csv(results_filename, delimiter=";", names=["language", "year", "estimator", "exponent", "fmin", "kernel size"])

	df = df.drop_duplicates()

	print(df)

	#estimator = "PDF Regression top 5000 types"
	estimator = "Better Try 30000 top words. Powerlaw beta fmin r_max on f(z)"

	df = df[(df['language'] == language_code) & (df['estimator'] == estimator)]

	print(df)

	df = df[df["year"] > min_year]

	print(df.head())

	plt.scatter(df["year"], df["exponent"], s=5)
	plt.title("Zipf exponent over time {}. Google 1grams\n{}".format(language_code, estimator).replace("_", " ").title())
	plt.savefig("images/freq_zipf_over_time_all_years_{}_{}_min_{}.png".format(language_code, estimator, min_year).replace(" ", "_").lower())
	plt.show()

def plot_all_years_fmin(min_year = 1800, language_code="eng-fiction"):

	results_filename = "results/powerlaw_zipf_results_by_year.csv"

	df = pd.read_csv(results_filename, delimiter=";", names=["language", "year", "estimator", "exponent", "fmin", "kernel size"])

	df = df.drop_duplicates()

	print(df)

	#estimator = "PDF Regression top 5000 types"
	estimator = "Better Try 30000 top words. Powerlaw beta fmin r_max on f(z)"

	df = df[(df['language'] == language_code) & (df['estimator'] == estimator)]

	print(df)

	df = df[df["year"] > min_year]

	print(df.head())

	plt.scatter(df["year"], df["fmin"], s=5)
	plt.title("Zipf exponent over time {}. Google 1grams\n{}".format(language_code, estimator).replace("_", " ").title())
	plt.savefig("images/fmin_zipf_over_time_all_years_{}_{}_min_{}.png".format(language_code, estimator, min_year).replace(" ", "_").lower())
	plt.show()


def plot_all_years_kernel_size(min_year = 1800, language_code="eng-fiction"):

	results_filename = "results/powerlaw_zipf_results_by_year.csv"

	df = pd.read_csv(results_filename, delimiter=";", names=["language", "year", "estimator", "exponent", "fmin", "kernel size"])

	df = df.drop_duplicates()

	print(df)

	#estimator = "PDF Regression top 5000 types"
	estimator = "Better Try 30000 top words. Powerlaw beta fmin r_max on f(z)"

	df = df[(df['language'] == language_code) & (df['estimator'] == estimator)]

	print(df)

	df = df[df["year"] > min_year]

	print(df.head())

	plt.scatter(df["year"], df["kernel size"], s=5)
	plt.title("Kernel Size over time {}. Google 1grams\n{}".format(language_code, estimator).replace("_", " ").title())
	plt.savefig("images/kernel_size_zipf_over_time_all_years_{}_{}_min_{}.png".format(language_code, estimator, min_year).replace(" ", "_").lower())
	plt.show()


def plot_all_languages():

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
		plot_all_years_kernel_size(min_year=1900, language_code=language_code) 

def calculate_theroetical_proportion_in_kernel_from_zipf(kernel_size, beta):

	alpha = 1/(beta-1)

	ps = get_probabilities_zeta_power_law(alpha, kernel_size)
	print(ps)
	total_p = sum(ps)
	print(total_p)




calculate_theroetical_proportion_in_kernel_from_zipf(1000, 1.8)