

import csv
import zipfile
import pandas as pd

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

			elif row_count < 100:
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


analyse_standardised_gutenberg()