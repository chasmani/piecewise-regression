import matplotlib.pyplot as plt
import pandas as pd

def plot_volume_vs_rank_frequency_one_year(year=2001):

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



	df.plot.scatter(x="rank", y="volume_count", s=1, alpha=0.2)

	plt.xscale("log")
	plt.yscale("log")

	plt.xlabel("frequency rank")
	plt.ylabel("volumes appeared in")
	plt.title("Frequency Rank agaist Volume Count of Google 1-grams {}".format(year))

	plt.savefig("images/volume_vs_freq_rank_{}.png".format(year))



plot_volume_vs_rank_frequency_one_year()