
import matplotlib.pyplot as plt
import pandas as pd

def plot_counts_dist_one_year(year=2001):

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

	df_c = df.groupby(by="match_count").count().reset_index()

	print(df_c.head())

	print(df_c.tail())

	

	df_c.plot.scatter(x="match_count", y="ngram", s=5)

	plt.xscale("log")
	plt.ylim(0.8, 1000000)

	plt.yscale("log")

	plt.xlabel("number of occurences of token, z")
	plt.ylabel("$\phi(z)$ - Number of tokens with this count")
	plt.title("Counts Distribution of Google 1-grams {}".format(year))

	plt.savefig("images/year_words_only_count_dist_plot_{}.png".format(year))


	


plot_counts_dist_one_year()