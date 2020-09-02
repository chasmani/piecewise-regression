install.packages("segmented")
library(segmented)
library(ggplot2)

# Set working directory if running from r studio
tryCatch({
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}, warning = function(w) {
}, error = function(e) {
}, finally = {
})

# Set working directoy if running from terminal
tryCatch({
  setwd(getSrcDirectory()[1])
}, warning = function(w) {
}, error = function(e) {
}, finally = {
})


analyseYear <- function(languageCode, year, outputFilename) {
  # Get Data

  inputFilename <- sprintf("ngrams_by_year/%s-%d.csv", languageCode, year)
  columns <- c("ngram", "year", "match_count", "volume_count", "language")
  myData <- read.csv(inputFilename, sep=";", col.names = columns) 
  myData <- subset(myData, match_count > 20)
  
  # Get ranks
  myData <- myData[order(myData$match_count, decreasing=TRUE),]
  myData$rank = 1:nrow(myData)
  
  # Log transforms
  myData$logrank <- log(myData$rank)
  myData$logfreq <- log(myData$match_count)
  
  # Linear model, needede for segmented analysis
  out.lm<-lm(logfreq~logrank,data=myData)
  
  # Breakpoint regression analysis with 2 breakpoints
  o2 <-segmented(out.lm, npsi=2, psi=c(9, 10.5)) 
  
  o2.breakpoint1 <- o2$psi["psi1.logrank", "Est."]
  o2.breakpoint2 <- o2$psi["psi2.logrank", "Est."]
  o2.alpha <- o2$coefficients["logrank"]
  o2.beta1 <- o2$coefficients["U1.logrank"]
  o2.beta2 <- o2$coefficients["U2.logrank"]
  o2.intercept <- o2$coefficients["(Intercept)"]
  
  print(o2$psi)
  print(o2$coefficients)
  
  
  csvData <- data.frame("Double Breakpoint R", languageCode, year, o2.breakpoint1, o2.breakpoint2, o2.alpha, o2.beta1, o2.beta2, o2.intercept)
  
  write.table(csvData,  
             file=outputFilename, 
             append = T, 
             sep=';', 
             row.names=F, 
             col.names=F )
  
}	

  
outputFilename <- "results/r_double_breakpoint_analysis_b.csv"
csvHeaders <- data.frame("Analysis", "Langauge Code", "Year", "Breakpoint 1", "Breakpoint 2", "alpha", "beta1", "beta2", "intercept")
write.table(csvHeaders, file=outputFilename, append = T, sep=';', row.names=F, col.names=F )

languageCode <- "eng-gb"
years <- 1995:2000

languages <- c("chi-sim", "eng-fiction",  "fre-all", "ger-all", "heb-all", "ita-all", "rus-all", "spa-all", "eng-us-all", "eng-all", "eng-gb")

languages <- c("eng-us-all", "eng-all", "eng-gb", "rus-all")


analyseYear("eng-us-all", 1968, outputFilename)

for (languageCode in languages){
  for (year in years){
    print(languageCode)
    print(year)
    
    # Trycatch because segmented sometimes doesn't work 
    tryCatch({
      analyseYear(languageCode, year, outputFilename)
    }, warning = function(w) {
    }, error = function(e) {
      print("Error")
      print(e)
      errorData <- data.frame("Double Breakpoint R", languageCode, year, "Error")
      write.table(errorData, file=outputFilename, append = T, sep=';', row.names=F, col.names=F )    
    }, finally = {
    })
  }
}




