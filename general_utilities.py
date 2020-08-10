import csv

def append_to_csv(csv_list, output_filename):
	with open(output_filename, "a", newline='') as fp:
		a = csv.writer(fp, delimiter=';')
		data=[csv_list]
		a.writerows(data)