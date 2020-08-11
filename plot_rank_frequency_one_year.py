
import matplotlib.pyplot as plt
import pandas as pd

def plot_rank_frequency_one_year(year=2001):

	df = pd.read_csv("data_by_year/{}.csv".format(year), delimiter=";", names=["ngram", "year", "match_count", "volume_count"])

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

	df = df.sort_values(by="match_count", ascending=False).reset_index(drop=True).reset_index()

	df["rank"] = df["index"] + 1

	print(df.head())



	df.plot.scatter(x="rank", y="match_count", s=5)

	plt.xscale("log")
	plt.yscale("log")

	plt.xlabel("rank")
	plt.ylabel("frequency")
	plt.title("Rank-frequency of Google 1-grams {}".format(year))

	plt.savefig("images/year_words_only_rank_freq_plot_{}.png".format(year))



plot_rank_frequency_one_year()