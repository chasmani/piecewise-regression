
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import pandas as pd

from general_utilities import append_to_csv

def generate_double_breakpoint_data():

	alpha = -4
	beta_1 = -2
	beta_2 = -2
	intercept = 100
	breakpoint_1 = 7
	breakpoint_2 = 10

	n_points = 200

	xx = np.linspace(0, 20, n_points)

	yy = intercept + alpha*xx + beta_1 * np.maximum(xx - breakpoint_1, 0) + beta_2 * np.maximum(xx - breakpoint_2, 0)  + np.random.normal(size=n_points)

	return xx, yy

def breakpoint_fit(xx, yy, current_breakpoint_1=8, current_breakpoint_2=10.5):

	A = xx
	U1 = (xx - current_breakpoint_1) * np.heaviside(xx- current_breakpoint_1, 1) 
	V1 = np.heaviside(xx - current_breakpoint_1, 1)

	U2 = (xx - current_breakpoint_2) * np.heaviside(xx- current_breakpoint_2, 1)
	V2 = np.heaviside(xx - current_breakpoint_2, 1)

	
	Z = np.array([A, U1 , U2, V1, V2]).T

	Z = sm.add_constant(Z)

	results = sm.OLS(endog=yy, exog=Z).fit()

	
	beta_1_hat = results.params[2]
	beta_2_hat = results.params[3]
	gamma_1_hat = results.params[4]
	gamma_2_hat = results.params[5]
	
	next_breakpoint_1 = current_breakpoint_1 - gamma_1_hat/beta_1_hat
	
	next_breakpoint_2 = current_breakpoint_2 - gamma_2_hat/beta_2_hat
	
	return next_breakpoint_1, next_breakpoint_2, results.params
	

def double_breakpoint_iterate():


	xx, yy = generate_double_breakpoint_data()
	current_breakpoint_1 = 9
	current_breakpoint_2 = 9.1




	for i in range(6):
		current_breakpoint_1, current_breakpoint_2, params = breakpoint_fit(xx, yy, current_breakpoint_1, current_breakpoint_2)
		print(current_breakpoint_1, current_breakpoint_2)

	intercept = params[0]
	alpha_hat = params[1]
	beta_1_hat = params[2]
	beta_2_hat = params[3]
	breakpoint_1_hat = current_breakpoint_1
	breakpoint_2_hat = current_breakpoint_2

	yy_hats = intercept + alpha_hat*xx + beta_1_hat * np.maximum(xx - breakpoint_1_hat, 0) + beta_2_hat * np.maximum(xx - breakpoint_2_hat, 0)

	plt.plot(xx, yy_hats, linewidth=2, color="purple", linestyle="dashed")
	plt.scatter(xx, yy, s=4, color="green")

	plt.show()
	
def double_breakpoint_analysis(xx, yy, current_breakpoint_1=9, current_breakpoint_2=10):

	delta = 1

	run_count = 0
	while delta > 0.001 and run_count<10:
		last_breakpoint_1, last_breakpoint_2 = current_breakpoint_1, current_breakpoint_2
		current_breakpoint_1, current_breakpoint_2, params = breakpoint_fit(xx, yy, current_breakpoint_1, current_breakpoint_2)
		delta = max(abs(current_breakpoint_1- last_breakpoint_1), abs(current_breakpoint_2 - last_breakpoint_2))
		run_count += 1


	return current_breakpoint_1, current_breakpoint_2, params

def double_breakpoint_test(year, language):

	input_filename = "ngrams_by_year/{}-{}.csv".format(language, year)
	df = pd.read_csv(input_filename, sep=";", names=["word", "year", "match_count", "volume_count", "language"])

	df = df[(df['match_count'] >20)]

	# Sort df
	df = df.sort_values(by="match_count", ascending=False)

	df = df.reset_index()

	df["rank"] = df.index + 1
	df["logrank"] = np.log(df["rank"])
	df["logcount"] = np.log(df["match_count"])

	xx = df["logrank"].to_numpy()
	yy = df["logcount"].to_numpy()

	breakpoint_1, breakpoint_2, params = double_breakpoint_analysis(xx, yy, 9, 10)

	intercept = params[0]
	alpha_hat = params[1]
	beta_1_hat = params[2]
	beta_2_hat = params[3]
	breakpoint_1_hat = breakpoint_1
	breakpoint_2_hat = breakpoint_2

	yy_hats = intercept + alpha_hat*xx + beta_1_hat * np.maximum(xx - breakpoint_1_hat, 0) + beta_2_hat * np.maximum(xx - breakpoint_2_hat, 0)

	plt.plot(xx, yy_hats, linewidth=2, color="red", linestyle="dashed")
	plt.axvline(breakpoint_1)
	plt.axvline(breakpoint_2)
	plt.scatter(xx, yy, s=2, color="lime")
	
	plt.xlabel("log(rank)")
	plt.ylabel("log(count)")

	plt.title("Double Breakpoint Analysis {} {}".format(year, language))

	plt.savefig("images/double_breakpoint_python_{}-{}.png".format(language, year))

	plt.show()

LANGAUGE_CODES = [
	"chi-sim", 
	"eng-fiction",  
	"fre-all", 
	"ger-all", 
	"heb-all", 
	"ita-all", 
	"rus-all", 
	"spa-all", 
	"eng-us-all", 
	"eng-all",
	"eng-gb"
]

def double_breakpoint_analysis_year_language(year, language):

	input_filename = "ngrams_by_year/{}-{}.csv".format(language, year)
	df = pd.read_csv(input_filename, sep=";", names=["word", "year", "match_count", "volume_count", "language"])

	df = df[(df['match_count'] >20)]

	# Sort df
	df = df.sort_values(by="match_count", ascending=False)

	df = df.reset_index()

	df["rank"] = df.index + 1
	df["logrank"] = np.log(df["rank"])
	df["logcount"] = np.log(df["match_count"])

	xx = df["logrank"].to_numpy()
	yy = df["logcount"].to_numpy()

	breakpoint_1, breakpoint_2, params = double_breakpoint_analysis(xx, yy, 9, 10)



	intercept = params[0]
	alpha_hat = params[1]
	beta_1_hat = params[2]
	beta_2_hat = params[3]
	breakpoint_1_hat = breakpoint_1
	breakpoint_2_hat = breakpoint_2

	print(alpha_hat, breakpoint_1, breakpoint_2, params)

	return [alpha_hat, breakpoint_1_hat, breakpoint_2_hat, intercept, beta_1_hat, beta_2_hat]


def breakpoint_analysis_on_all():

	output_filename = "results/double_breakpoint_python_historical.csv"

	for language in LANGAUGE_CODES:
		for year in range(1800, 2000):
			try:
				print("Working on {}-{}".format(language, year))
				results = double_breakpoint_analysis_year_language(year, language)
				csv_row = ["Breakpoint Analysis Python", year, language] + results
				append_to_csv(csv_row, output_filename)
			except Exception as e:
				print(str(e))
				csv_row = ["Breakpoint Analysis Python", year, language, "Error", str(e)]





if __name__=="__main__":
	breakpoint_analysis_on_all()