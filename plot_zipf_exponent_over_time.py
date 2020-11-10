
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_all_years_clauset(min_year = 1800, language_code="eng-fiction"):

	results_filename = "results/zipf_results_by_year.csv"

	df = pd.read_csv(results_filename, delimiter=";", names=["language", "year", "estimator", "exponent"])

	df = df.drop_duplicates()

	#estimator = "PDF Regression top 5000 types"
	estimator = "Clauset top 5000 types"

	df = df[(df['language'] == language_code) & (df['estimator'] == estimator)]

	df = df[df["year"] > min_year]

	print(df.head())

	plt.scatter(df["year"], df["exponent"], s=5)
	plt.title("Zipf exponent over time {}. Google 1grams\n{}".format(language_code, estimator).replace("_", " ").title())
	plt.savefig("images/zipf_over_time_all_years_{}_{}_min_{}.png".format(language_code, estimator, min_year).replace(" ", "_").lower())
	plt.show()


def plot_kernel_size_over_time_breakpoints(min_year=1800, language_code="fre-all"):

	results_filename = "results/r_double_breakpoint_analysis_b.csv"

	df = pd.read_csv(results_filename, delimiter=";", header=0)

	df = df[(df['Langauge Code'] == language_code) & (df['Analysis'] == "Double Breakpoint R") & (df["Breakpoint 1"] != "Error")]

	df["Year"] = pd.to_numeric(df["Year"])
	df["Breakpoint 1"] = pd.to_numeric(df["Breakpoint 1"])
	df["Kernel Size"] = np.exp(df["Breakpoint 1"])

	df = df[df["Year"] > min_year]
	df.drop_duplicates(subset="Year")

	df.plot.scatter(x="Year", y="Kernel Size", s=3)

	plt.title("Kernel Size Over Time {}. \nBreakpoint Analysis on Google 1grams".format(language_code))
	plt.savefig("images/breakpoint_kernel_size_{}_min_{}.png".format(language_code, min_year).replace(" ", "_").lower())
	plt.show()


def plot_zipf_exponent_over_time_breakpoints(min_year=1800, language_code="fre-all"):

	results_filename = "results/r_double_breakpoint_analysis_b.csv"

	df = pd.read_csv(results_filename, delimiter=";", header=0)

	df = df[(df['Langauge Code'] == language_code) & (df['Analysis'] == "Double Breakpoint R") & (df["Breakpoint 1"] != "Error")]



	print(df["Year"])

	print(df.dtypes)

	df["Year"] = pd.to_numeric(df["Year"])
	df["Breakpoint 1"] = pd.to_numeric(df["Breakpoint 1"])
	df["Kernel Exponent"] = -1*pd.to_numeric(df["alpha"])

	

	df = df[df["Year"] > min_year]

	df.drop_duplicates(subset="Year")

	df.plot.scatter(x="Year", y="Kernel Exponent", s=3)

	plt.title("Kernel Zipf Exponent Over Time {}. \nBreakpoint Analysis on Google 1grams".format(language_code))
	plt.savefig("images/breakpoint_zipf_exponent_{}_min_{}.png".format(language_code, min_year).replace(" ", "_").lower())
	plt.show()


def plot_kernel_size_and_exponent_over_time_breakpoints(min_year=1800, language_code="fre-all"):


	if language_code == "eng-gb":
		results_filename = "results/r_double_breakpoint_analysis.csv"		
	else:
		results_filename = "results/r_double_breakpoint_analysis_b.csv"

	df = pd.read_csv(results_filename, delimiter=";", header=0)

	df = df[(df['Langauge Code'] == language_code) & (df['Analysis'] == "Double Breakpoint R") & (df["Breakpoint 1"] != "Error")]



	print(df["Year"])

	print(df.dtypes)

	df["Year"] = pd.to_numeric(df["Year"])
	df["Breakpoint 1"] = pd.to_numeric(df["Breakpoint 1"])
	df["Kernel Exponent"] = -1*pd.to_numeric(df["alpha"])
	df["Kernel Size"] = np.exp(df["Breakpoint 1"])

	df = df[df["Year"] > min_year]
	if language_code == "eng-us-all":
		df = df[df["Kernel Exponent"] < 1.05]
		#df = df[df["Kernel Size"] > min_year]
		



	df.drop_duplicates(subset="Year")


	f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

	df.plot.scatter(x="Year", y="Kernel Exponent", s=3, ax = ax1)
	df.plot.scatter(x="Year", y="Kernel Size", s=3, ax = ax2)
	
	"""
	ax1.axvspan(1914, 1918, alpha=0.5, color='red')
	ax1.axvspan(1939, 1945, alpha=0.5, color='red')
	ax2.axvspan(1914, 1918, alpha=0.5, color='red')
	ax2.axvspan(1939, 1945, alpha=0.5, color='red')
	"""


	plt.suptitle("Kernel Over Time {}. \nBreakpoint Analysis on Google 1grams".format(language_code))
	plt.savefig("images/breakpoint_kernel_both_{}_min_{}.png".format(language_code, min_year).replace(" ", "_").lower())
	plt.show()


def plot_kernel_size_vs_exponent(min_year=1800, language_code="fre-all"):


	if language_code == "eng-gb":
		results_filename = "results/r_double_breakpoint_analysis.csv"		
	else:
		results_filename = "results/r_double_breakpoint_analysis_b.csv"

	df = pd.read_csv(results_filename, delimiter=";", header=0)

	df = df[(df['Langauge Code'] == language_code) & (df['Analysis'] == "Double Breakpoint R") & (df["Breakpoint 1"] != "Error")]



	print(df["Year"])

	print(df.dtypes)

	df["Year"] = pd.to_numeric(df["Year"])
	df["Breakpoint 1"] = pd.to_numeric(df["Breakpoint 1"])
	df["Kernel Exponent"] = -1*pd.to_numeric(df["alpha"])
	df["Kernel Size"] = np.exp(df["Breakpoint 1"])

	df = df[df["Year"] > min_year]
	if language_code == "eng-us-all":
		df = df[df["Kernel Exponent"] < 1.05]
		#df = df[df["Kernel Size"] > min_year]
		



	df.drop_duplicates(subset="Year")

	df.plot.scatter(x="Kernel Size", y="Kernel Exponent", s=3, c="Year")
	
	"""
	ax1.axvspan(1914, 1918, alpha=0.5, color='red')
	ax1.axvspan(1939, 1945, alpha=0.5, color='red')
	ax2.axvspan(1914, 1918, alpha=0.5, color='red')
	ax2.axvspan(1939, 1945, alpha=0.5, color='red')
	"""


	plt.suptitle("Kernel Size vs Exponent {}. \nBreakpoint Analysis on Google 1grams".format(language_code))
	plt.savefig("images/kernel_size_vs_exponent_{}_min_{}.png".format(language_code, min_year).replace(" ", "_").lower())
	plt.show()




FREANCH_WARS = [
	(1803, 1815),
	(1914, 1918),
	(1939, 1945)
]

GERMAN_WARS = [
	(1805, 1815),
	(1870, 1871),
	(1914, 1918),
	(1939, 1945)
]


BRITISH_WARS = [
	(1914, 1918),
	(1939, 1945)
]








def plot_all_languages():




	for language_code in [
		"eng-gb",
		"eng-fiction",
		"eng-all", 
		"eng-us-all", 
		"fre-all", 
		"ger-all", 
		"heb-all", 
		"ita-all", 
		"rus-all", 
		"spa-all"]:
		plot_kernel_size_vs_exponent(min_year=1800, language_code=language_code) 

plot_kernel_size_vs_exponent(min_year=1800, language_code="eng-fiction")