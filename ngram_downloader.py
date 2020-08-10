

import requests
import os
#
URL = "http://storage.googleapis.com/books/ngrams/books/googlebooks-{}-all-1gram-20120701-{}.gz"

LANGUAGES = [
	"eng",
	"eng-us",
	"eng-gb",
	"eng-fiction",
	"chi-sim",
	"fre",
	"ger",
	"heb",
	"ita",
	"rus",
	"spa"
]
FILE_AFFIXES = ["0","1","2","3","4","5","6","7","8","9","a","b","c","d","e","f","g","h","i","j","k","l",
	"m","n","o","other","p","pos","punctuation","q","r","s","t","u","v","w","x","y","z"]

def download_all():

	for language in LANGUAGES:
		for affix in FILE_AFFIXES:
			url = URL.format(language, affix)
			filename = "ngram_raw_data/" + url.split("/")[-1]
			if os.path.isfile(filename):
				print("Already got: \t", filename)
			else:
				print("Downloading: \t", filename)
				r = requests.get(url, allow_redirects=True)
				open(filename, 'wb').write(r.content)						



download_all()