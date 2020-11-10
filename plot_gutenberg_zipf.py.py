

import csv
import zipfile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import zipfanalysis

from general_utilities import append_to_csv


def analyse_standardised_gutenberg():

	metadata_file = "gutenberg_standard/SPGC-metadata-2018-07-18.csv"
	output_file = "results/gutenberg_analysis.csv"

	with open(metadata_file, newline='') as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		row_count = 0
		for row in reader:

			print("Working row {} of 57714 . . . ".format(row_count))
			
			if row_count == 0:
				csv_row = row + ["alpha_pdf", "total_words", "unique_words"]
				append_to_csv(csv_row, output_file)

			elif row_count > 0:
				gut_id = row[0]
				title = row[1]
				author = row[2]
				author_birth_year = row[3]
				author_birth_death = row[4]
				language = row[5]
				downloads = row[6]
				subjects = row[7]
				filetype = row[8]

				try:

					zf = zipfile.ZipFile('gutenberg_standard/SPGC-counts-2018-07-18.zip') 
					df = pd.read_csv(zf.open('SPGC-counts-2018-07-18/{}_counts.txt'.format(gut_id)), delimiter="\t", names=["ngram", "count"])
					n = df["count"].to_numpy()

					alpha = zipfanalysis.estimators.ols_regression_pdf.ols_regression_pdf_estimator(n)
					total_words = df["count"].sum()
					unique_words = len(n)

					csv_row = row + [alpha, total_words, unique_words]

				except Exception as e:
					csv_row = row + [str(e)]

				print(csv_row)
				append_to_csv(csv_row, output_file)


			row_count += 1


def plot_gutenberg_zipf():

	min_year = 1800

	input_file = "results/gutenberg_analysis.csv"

	df = pd.read_csv(input_file, sep=";", names=["gut_id", "title", "author", "author_birth_year", 
			"author_birth_death", "language", "downloads", "subjects", "filetype", "alpha", "total_words", "unique_words"])

	df = df[(df['language'] == "['en']")]

	df["author_birth_year"] = pd.to_numeric(df["author_birth_year"], errors="coerce")
	df["alpha"] = pd.to_numeric(df["alpha"], errors="coerce")




	df = df[(~df['author_birth_year'].isna())]
	df = df[(~df['alpha'].isna())]

	df = df[(df['author_birth_year'] > min_year)]

	df['Author Birth Date'] = pd.cut(df['author_birth_year'], bins=range(min_year, 2000, 10))
	sns.boxplot(x="Author Birth Date", y="alpha", data=df)

	plt.xticks(rotation=45, ha="right")
	

	plt.title("Author Birth Date and Zipf Exponent of Books\nProject Gutenberg")
	plt.tight_layout()

	plt.savefig("images/gutenberg_author_birth_date_alpha_boxplots_min_year_{}.png".format(min_year))

	plt.show()

def plot_gutenberg_zipf_grouped_by_length():

	min_year = 1800

	input_file = "results/gutenberg_analysis.csv"

	df = pd.read_csv(input_file, sep=";", names=["gut_id", "title", "author", "author_birth_year", 
			"author_birth_death", "language", "downloads", "subjects", "filetype", "alpha", "total_words", "unique_words"])

	df = df[(df['language'] == "['en']")]

	df["author_birth_year"] = pd.to_numeric(df["author_birth_year"], errors="coerce")
	df["alpha"] = pd.to_numeric(df["alpha"], errors="coerce")
	df["unique_words"] = pd.to_numeric(df["unique_words"], errors="coerce")




	df = df[(~df['author_birth_year'].isna())]
	df = df[(~df['alpha'].isna())]
	df = df[(~df['unique_words'].isna())]

	df = df[(df['author_birth_year'] > min_year)]



	df['Author Birth Date'] = pd.cut(df['author_birth_year'], bins=range(min_year, 2000, 10))

	df["Unique Words Bins"] = pd.cut(df['unique_words'], bins=range(0, 10000, 1000))

	df["Author Birth Date Mean"] = df["Author Birth Date"].apply(lambda x: x.mid)

	df_grouped = df.groupby(["Author Birth Date Mean", "Unique Words Bins"]).mean().reset_index()

	print(df_grouped)

	sns.lineplot(x="Author Birth Date Mean", y="alpha", hue="Unique Words Bins", data=df_grouped)

	plt.xticks(rotation=45, ha="right")

	
	plt.title("Author Birth Date and Zipf Exponent of Books with Unique Word Bins\nProject Gutenberg")
	#plt.tight_layout()

	plt.savefig("images/gutenberg_author_birth_date_alpha_boxplots_unique_word_bins_min_year_{}.png".format(min_year))

	plt.show()
	

plot_gutenberg_zipf_grouped_by_length()