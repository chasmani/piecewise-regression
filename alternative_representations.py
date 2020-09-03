import pandas as pd
import matplotlib.pyplot as plt

def see_all_representations(language="eng-gb", year=1887):


	input_filename = "ngrams_by_year/{}-{}.csv".format(language, year)
	df = pd.read_csv(input_filename, sep=";", names=["word", "year", "match_count", "volume_count", "language"])
	df = df.sort_values(by="match_count", ascending=False)
	df = df.reset_index()
	df["rank"] = df.index + 1

	df = df[(df['match_count'] >20)]


	"""
	# Frequency rank representation
	df.plot.scatter(x="rank", y="match_count", label="Frequency-Rank", s=2)
	plt.xscale("log")
	plt.yscale("log")

	plt.show()
	"""

	# CCDF of frequency rank
	# Version 1 - one point for every match
	# Version 2 - one point for every rank
	# Version 3 - one point for every tied rank
	
	# Version 1 is very computationally expensive for large datasets

	# Version 2
	df["ccdf_match_count"] = df["match_count"].sum() - df["match_count"].cumsum().shift(1)
	#df["ccdf_match_count"][0] = 1

	#/df["match_count"].sum()


	df.plot.scatter(x="rank", y="ccdf_match_count", s=2)
	plt.xscale("log")
	plt.yscale("log")

	plt.ylim(df["ccdf_match_count"].min()*0.5, df["ccdf_match_count"].max()*1.5)

	plt.show()
	print(df)


see_all_representations()
