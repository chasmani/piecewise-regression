import os
import csv
import gzip

from general_utilities import append_to_csv

def extract_language_files_to_csvs_by_year(language_code="eng-gb"):
	zipped_files = os.listdir("ngram_raw_data")
	this_language_zipped_files = [f for f in zipped_files if language_code in f]
	this_language_zipped_files_ready = [f for f in this_language_zipped_files if "punctuation" not in f]

	for zip_filename in sorted(this_language_zipped_files_ready):
		print("Working . . . ", zip_filename)

		with gzip.open("ngram_raw_data/"+zip_filename, mode="rt") as f:
			
			reader = csv.reader(f, delimiter="\t")
			
			for row in reader:
				ngram = row[0]
				year=row[1]
				match_count = row[2]
				volume_count=row[3]
				row.append(language_code)
				# Don't include words tagged with _ADJ, _NOUN etc 
				if "_" not in ngram:
					append_to_csv(row, "ngrams_by_year/{}-{}.csv".format(language_code, year))
			

if __name__=="__main__":
	#for language_code in ["chi-sim", "eng-all", "eng-us-all", "fre-all", "ger-all", "heb-all", "ita-all", "rus-all", "spa-all"]:
	for language_code in ["ita-all", "rus-all", "spa-all"]:	
		extract_language_files_to_csvs_by_year(language_code)