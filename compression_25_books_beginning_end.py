

import zipfanalysis

import os
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt

import re
import unidecode
from collections import Counter
from general_utilities import append_to_csv

import lzma

def compress_books():

	directory = "gutenberg_manual"

	# 1. Get all files
	files = [f for f in listdir(directory)]
	print(files)
	for file in files:
		print(file)
		file_dir = "gutenberg_manual/" + file
		bookname = file.replace("_", " ").title()
		with open (file_dir, 'r', errors="ignore") as f:
			content = f.read()
			
			total_words = len(content)
			print(total_words)
			first_half = content[:int(total_words/2)]
			second_half = content[int(total_words/2):]

			first_half_bytes = str.encode(first_half)
			obj_first = lzma.LZMAFile(file_dir + "_half_1.xz", mode="wb")
			obj_first.write(first_half_bytes)
			obj_first.close()	

			second_half_bytes = str.encode(second_half)
			obj = lzma.LZMAFile(file_dir + "_half_2.xz", mode="wb")
			obj.write(second_half_bytes)
			obj.close()	




			"""

			data = b'Welcome to TutorialsPoint'
>>>obj = lzma.LZMAFile("test.xz", mode="wb")
>>>obj.write(data)
>>>obj.close()
"""






	# 2. Remove project gutenberg bits

if __name__=="__main__":
	compress_books()