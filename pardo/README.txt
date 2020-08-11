This folder contains the data, scripts and the final report of my project. It consists of four parts:

*OriginalData: The folder where you can find the data used for the project and the scripts used to get it from Google's available data. It has two subfolders:
	
	*Scripts: It includes two scripts:
		*ngrams.sh: Script found on the internet to download the n-grams 					authomatically.
		*ZipfScript: Script to create the final tables that were used.

	*Tables: It includes the six tables of word frequencies used in the project, one 		for each language.
		Their structure is as follows:
			
			word 1800 1801 . . . 2011 2012 total
			a     101 150  . . . 1000 1153 10570
			aa     0   6   . . .  100  106  560
			.      .   .   . . .   .    .    .
			.      .   .   . . .   .    .    .
			.      .   .   . . .   .    .    .
		
		Where total is the sum of the frequencies for all years.
	
*Analysis: Contains the R scripts with which I analyzed and plotted the data. There are four of them:

	*Main.R: The main file in which all the code is included. Although it is 			documented, it is not very organized and it can be difficult to make sense 		of some of its sections (especially when data created somewhere else is 		imported). For this reason I created two independent scripts to make the 		data analysis more clear. However the code for plotting the final results		 is found only in this script. Almost all the code included in this script 		was used to obtain the results shown in the report.

	*SegmentedAndPowerLaw.R: A script to obtain the results when a bootstrap method is 				apply to the data obtains through both the power law fit 				and segmented regression methods. The data is saved as 					an .Rda file, making it available to be imported in Main.R 				to generate the plots.
	
	*BootstrapKernel.R: A script where the results for the fixed kernel with bootstrap 				are obtained. The data is saved as 							an .Rda file, making it available to be imported in Main.R 				to generate the plots.

	*ZipfYears.R: The initial file for the data analysis of the project. Really 				unorganized and scarcely documented. However, it includes a lot of 			experiments that were not included in the final report and some 			that I mentioned (seeing the data as time series and looking for 			correlations). Some useful code snippets might be found in this 			file.

*FinalData: Includes the data generated using the scripts mentioned above (.Rda 			files) which can be imported into R to generate the plots or to do 			further investigations. The files whose names start with Boots 				were created using BotstrapKernel.R, while the ones starting with 			coef were created by SegmentedAndPowerLaw.R.

*Report: It contains the necessary files to generate the the report (.tex, .bib, plots and 	a file with the Latex template I used (revcoles.cls)).



	

