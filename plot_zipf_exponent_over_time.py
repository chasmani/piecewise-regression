
import pandas as pd
import matplotlib.pyplot as plt

def plot_all_years(min_year = 1800):

	results_filename = "results/zipf_results_by_year.csv"

	df = pd.read_csv(results_filename, delimiter=";", names=["language", "year", "estimator", "exponent"])

	df = df.drop_duplicates()

	language_code = "eng-fiction"
	estimator = "PDF Regression top 5000 types"
	estimator = "Clauset top 5000 types"

	df = df[(df['language'] == language_code) & (df['estimator'] == estimator)]

	df = df[df["year"] > min_year]

	print(df.head())

	plt.scatter(df["year"], df["exponent"], s=5)
	plt.title("Zipf exponent over time {}. Google 1grams\n{}".format(language_code, estimator).replace("_", " ").title())
	plt.savefig("images/zipf_over_time_all_years_{}_{}_min_{}.png".format(language_code, estimator, min_year).replace(" ", "_").lower())
	plt.show()

plot_all_years(1800)