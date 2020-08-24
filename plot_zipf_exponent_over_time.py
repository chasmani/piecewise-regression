
import pandas as pd
import matplotlib.pyplot as plt

def plot_all_years(min_year = 1800, language_code="eng-fiction"):

	results_filename = "results/zipf_results_by_year.csv"

	df = pd.read_csv(results_filename, delimiter=";", names=["language", "year", "estimator", "exponent"])

	df = df.drop_duplicates()

	#estimator = "PDF Regression top 5000 types"
	estimator = "Clauset top 5000 types"

	df = df[(df['language'] == language_code) & (df['estimator'] == estimator)]

	df = df[df["year"] > min_year]

	print(df.head())

	plt.scatter(df["year"], df["exponent"], s=5)
	plt.title("Zipf exponent over time {}. Google 1grams\n{}".format(language_code, estimator).replace("_", " ").title())
	plt.savefig("images/zipf_over_time_all_years_{}_{}_min_{}.png".format(language_code, estimator, min_year).replace(" ", "_").lower())
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
		plot_all_years(min_year=1900, language_code=language_code) 



plot_all_languages()