
import pandas as pd
import matplotlib.pyplot as plt

from general_utilities import append_to_csv

LANGAUGE_CODES = [
	"chi-sim", 
	"eng-fiction",  
	"fre-all", 
	"ger-all", 
	"heb-all", 
	"ita-all", 
	"rus-all", 
	"spa-all", 
	"eng-us-all", 
	"eng-all",
	"eng-gb"
]

def extract_mean_unique_words():

	output_filename = "results/mean_unique_words.csv"

	for language in LANGAUGE_CODES:
		for year in range(1800, 2000):
			input_filename = "ngrams_by_year/{}-{}.csv".format(language, year)
			df = pd.read_csv(input_filename, sep=";", names=["word", "year", "match_count", "volume_count", "language"])

			max_volume = df["volume_count"].max()
			total_unique_words = df["volume_count"].sum()
			mean_unique_words = total_unique_words/max_volume
			print(language, year, max_volume, total_unique_words, mean_unique_words)
			csv_row = ["Max Volume Count", language, year, max_volume, total_unique_words, mean_unique_words]


			append_to_csv(csv_row, output_filename)


def plot_mean_unique_words():

	input_filename = "results/mean_unique_words.csv"
	df = pd.read_csv(input_filename, sep=";", names=["analysis", "language", "year", "max_volume", "total_unique_words", "mean_unique_words"])

	print(df)

	for language in LANGAUGE_CODES:
		df_l = df[(df['language'] == language)]

		df_l.plot.scatter(x="year", y="mean_unique_words", s=3)
		plt.title("Mean Unique Words Per Year Ngrams {}".format(language))
		plt.savefig("images/mean_unique_words_{}.png".format(language))
		plt.show()




if __name__=="__main__":
	plot_mean_unique_words()