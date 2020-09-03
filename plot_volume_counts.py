
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

def extract_volumes_counts():

	output_filename = "results/max_volumes.csv"

	for language in LANGAUGE_CODES:
		for year in range(1800, 2000):
			input_filename = "ngrams_by_year/{}-{}.csv".format(language, year)
			df = pd.read_csv(input_filename, sep=";", names=["word", "year", "match_count", "volume_count", "language"])

			max_volume = df["volume_count"].max()
			print(language, year, max_volume)
			csv_row = ["Max Volume Count", language, year, max_volume]


			append_to_csv(csv_row, output_filename)

def plot_volumes():

	input_filename = "results/max_volumes.csv"
	df = pd.read_csv(input_filename, sep=";", names=["analysis", "language", "year", "max_volume"])

	print(df)

	for language in LANGAUGE_CODES:
		df_l = df[(df['language'] == language)]

		df_l.plot.scatter(x="year", y="max_volume", s=3)
		plt.title("Volumes Per Year Ngrams {}".format(language))
		plt.savefig("images/max_volumes_{}.png".format(language))
		plt.show()



if __name__=="__main__":
	plot_volumes()