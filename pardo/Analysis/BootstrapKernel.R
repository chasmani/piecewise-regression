library(igraph)
library(MASS)
library(stats)
library(SOAR)
library(segmented)
library(psych)
library(ggplot2)
library(poweRlaw)
library(reshape2)
library(gtools)

#Loads the data

#data_ger = read.csv("all_words_per_year_ger.txt")
#data_spa = read.csv("all_words_per_year_spa.txt")
#data_fre = read.csv("all_words_per_year_fre.txt")
data_ita = read.csv("all_words_per_year_ita.txt")
#data_eng = read.csv("all_words_per_year_eng.txt", quote = "", row.names = NULL, stringsAsFactors = FALSE)
#data_rus = read.csv("all_words_per_year_rus.txt", quote = "", row.names = NULL, stringsAsFactors = FALSE)


#Perfoms a bootstrap of the data. It selects 1000000 words (repetitions of items are accepted) with probabilities 
#corresponding to their frequencies in the original dataset.

bootstrapData<-function(data)
{
  n<-1000000

  data<-data[data!=0]
  data<-data[!is.na(data)]
  #Size of the array
  numTokens<-length(data)
  #total number of words
  totalItems<-sum(as.numeric(data))
  
  #The cumulative sums of the number of occurrences of words
  probs<-cumsum(as.numeric(data))
  #Random numbers to select words
  ranNums<-floor(runif(n,min=1,max=totalItems))
  #A data frame with the indexes of the words selected at random and the number of times they were selected
  indices<-as.data.frame(table(findInterval(ranNums,probs,all.inside=TRUE)))#sapply(ranNums,findPosition,cumVector=probs)))
  indices[[1]]<-as.numeric(as.character(indices[[1]]))
  newData<-rep(0,length(data))
  newData[indices[[1]]]<-indices[[2]]
  return(newData)
  
  
}

#Calculates Zipf's coefficient of the 5000 highest frequency words. It is done via a 
#linear regresion because we want to check that the evolution observed in the kernel lexicon
#is not due to size, and using Clauset's method the number of words taken into account would 
#change for each bootstrapped sample. Less words can be considered (e.g. 2000) to discard any 
#possible influence of the tail.
calculateCoefFix<-function(data)
{
  #Total number of words
  totalData<-sum(as.numeric(data))
  #Maximum frequency
  maxDoc<-max(as.numeric(data))
  #Eliminates data with 0 frequency
  data<-data[data!=0]
  #normalizes the data
  data<-data/totalData
  #Counts the number of items (number of different words)
  numWords<-length(data)
  #Creates a histogram out of the frequencies. (How many words are found in each interval of frequencies)
  data_hist<-try(hist(log(data),breaks = 30, plot=FALSE))
  #Captures an exception if the histogram cannot be created
  if(class(data_hist) == "try-error") {return(NA))}
  #Sorts the data
  x<-sort(data_hist$counts[data_hist$counts!=0],decreasing=FALSE)
  y<-sort(exp(data_hist$mids[data_hist$counts!=0]),decreasing=TRUE)
  
  x<-x[1:(length(x)-2)]
  y<-y[1:(length(y)-2)]
  
  #Logarithmic scale, where the linear regression will be performed
  x<-log(x)
  y<-log(y)
  #Linear regression
  if(length(data)!=0)
  {  
    m <- lm(y~x)
    value<-coef(m)[[2]]
  }
  else
    value<-NA
  return(-value)
}

#Selects the 5000 highest frequncy words.
fix<-function(data)
{
  data<-sort(data,decreasing=TRUE)
  return(data[1:5000])
}

#Generates 100 bootstrapped samples and calculates Zipf's coefficients for them.
statBootstrapFix<-function(data)
{
  newData<-replicate(100,bootstrapData(data))
  newData<-apply(newData,2,fix)
  totalStats<-apply(newData,2,calculateCoefFix)
  #Calculate the average of the bootstrap statistics
  #averages<-apply(totalStats,1,mean)
  #print("aaaaaa")
  return(totalStats)
}

#BootsEngFix<-apply(data_eng[2:209],2,statBootstrapFix)
#save(BootsEngFix.file="BootsEngFix.Rda")
#BootsGerFix<-apply(data_ger[2:209],2,statBootstrapFix)
#BootsRusFix<-apply(data_rus[2:209],2,statBootstrapFix)
#save(BootsRusFix.file="BootsEngRus.Rda")
#BootsSpaFix<-apply(data_spa[2:209],2,statBootstrapFix)
#BootsFreFix<-apply(data_fre[2:209],2,statBootstrapFix)

#It runs the bootstrap function for a particular data and saves the results in a file.
BootsItaFix<-apply(data_ita[4:209],2,statBootstrapFix)
save(BootsItaFix.file="BootsIta.Rda")