#!/usr/bin/env Rscript
#install.packages("igraph")
#install.packages("MASS")
#install.packages("stats")
#install.packages("SOAR")
#install.packages("segmented")
#install.packages("psych")
#install.packages("poweRlaw")
library(igraph)
library(MASS)
library(stats)
library(SOAR)
library(segmented)
library(psych)
library(ggplot2)
library(poweRlaw)

#data_ger = read.csv("all_words_per_year_ger.txt")
#data_spa = read.csv("all_words_per_year_spa.txt")
#data_fre = read.csv("all_words_per_year_fre.txt")
#data_ita = read.csv("all_words_per_year_ita.txt")
data_eng = read.csv("all_words_per_year_eng.txt", quote = "", row.names = NULL, stringsAsFactors = FALSE)
#data_rus = read.csv("all_words_per_year_rus.txt", quote = "", row.names = NULL, stringsAsFactors = FALSE)

#Determines whether the input 'var' represents a percentage
isPercentage<-function(var,total){
  if(var<2){
    var<-floor(total*var)
  }
  return(var)
}

#Calculates Zipf's coefficient of the kernel lexicon using Clauset's power law fitting method.

calculateCoefPow<-function(freqs,init,end)
{
  #Possible NA values are filtered out
  occur_per_year<-freqs[!is.na(freqs)]
  #years with no occcurrence of a word are omitted, since the logarithm will not be defined.
  occur_per_year<-occur_per_year[occur_per_year!=0]
  #The words are ranked
  occur_per_year<-sort(occur_per_year, decreasing = TRUE)
  limits<-list(init,end)
  limits<-lapply(limits,isPercentage,total=length(occur_per_year))
  occur_per_year<-occur_per_year[limits[[1]]:limits[[2]]]
  #Possible NA values are filtered out
  occur_per_year<-occur_per_year[!is.na(occur_per_year)]
  
  #Normalization
  total<-sum(occur_per_year)
  occur_per_yearNorm<-occur_per_year/total
  #The sequence of ranks is created
  x <- seq_along(occur_per_year)
  xmin<-0
  #The coefficients of the linear regressions are calculated
  if(length(occur_per_year)>xmin)
  {  
    m <- power.law.fit(occur_per_year,implementation='plfit')
    #Normalized
    m1<- power.law.fit(occur_per_yearNorm,implementation='plfit')
    value<-m$alpha
    value1<-m1$alpha
    test<-m$KS.stat
    test2<-m$KS.p
    xmin<-m$xmin
  }
  else
  {
    value<-NA
    test<-NA
    test2<-NA
    xmin<-NA
  }
  return(c(-value,-value1,total,test,test2,xmin))
  
}

data<-list(data_ger,data_spa,data_fre,data_ita,data_eng,data_rus)#,data_chi)
#results<-lapply(data,function(name){return(apply(name[2:209],2,calculateCoefPow,init=2,end=nrow(name)-1))})

#Calculates Zipf's coefficients of both the kernel and unlimited lexicons using a 
#segmented regression method. In this case it is not checked that the data actually fits a power law.
segmentedRegression<-function(data)
{
  #Total number of words
  totalData<-sum(data)
  #Highest frequency
  maxDoc<-max(data)
  #Eliminates data with 0 frequency
  data<-data[data!=0]
  #normalizes the data
  data<-data/totalData
  #Counts the number of items (number of different words)
  numWords<-length(data)
  #Creates a histogram out of the frequencies. (How many words are found in each interval of frequencies)
  data_hist<-try(hist(log(data),breaks = 20, plot=FALSE))
  if(class(data_hist) == "try-error") {return(c(NA,NA,NA,NA,NA))}
 
  x<-sort(data_hist$counts[data_hist$counts!=0],decreasing=FALSE)
  y<-sort(exp(data_hist$mids[data_hist$counts!=0]),decreasing=TRUE)
  x<-x[1:(length(x)-2)]
  y<-y[1:(length(y)-2)]
  #Where the power law stops holding, so is a good approximation for the discontinuity in 
  #the linear regression.
  aux<-power.law.fit(y,implementation="plfit")$xmin
  
  x<-log(x)
  y<-log(y)
  #The segmented regression is performed
  lin.mod<-lm(y~x)
  psi.init<-(log(aux)-lin.mod$coefficients[[1]])/lin.mod$coefficients[[2]]
  segmented.mod <- try(segmented(lin.mod, seg.Z = ~x, psi=7))
  if(class(segmented.mod) == "try-error") {
    #segmented.mod <-try(segmented(lin.mod, seg.Z = ~x, psi=7))
    #if(class(segmented.mod) == "try-error")
    return(c(NA,NA,NA,NA,NA))
  }
  #recovering the break point
  #Number of words in each column and its respective frequency
  a<-data_hist$counts[data_hist$counts>exp(segmented.mod$psi[2])]
  b<-data_hist$mids[data_hist$counts>exp(segmented.mod$psi[2])]
  b<-exp(b)
  b<-b*totalData
  #Total number of words in the kernel lexicon 
  c<-a*b
  coef1<-slope(segmented.mod)$x[1,1]
  coef2<-slope(segmented.mod)$x[2,1]
  psi<-segmented.mod$psi[2]
  return(c(sum(c)/totalData,-coef1,-coef2,aux,psi,psi.init))
}

bootstrapData<-function(data)
{
    n<-1000000
    data<-data[data!=0]
    data<-data[!is.na(data)]
    #Size of the array
    numTokens<-length(data)
    #total number of words
    totalItems<-sum(data)
    #The cumulative sums of the number of occurrences of words
    probs<-cumsum(data)
    #Random numbers to select words
    ranNums<-floor(runif(n,min=1,max=totalItems))
    #A data frame with the indexes of the words selected at random and the number of times they were selected
    indices<-as.data.frame(table(findInterval(ranNums,probs,all.inside=TRUE)))#sapply(ranNums,findPosition,cumVector=probs)))
    indices[[1]]<-as.numeric(as.character(indices[[1]]))
    newData<-rep(0,length(data))
    newData[indices[[1]]]<-indices[[2]]
    return(newData)
    
}

#Generates 100 bootstrapped samples and calculates Zipf's coefficients for them.
#Both methods (Power law fitting and segmented regression) can be used by uncommenting or commenting
# the second and third lines as desired (only one should be uncommented).
statBootstrap<-function(data)
{
  newData<-replicate(100,bootstrapData(data))
  totalStats<-apply(newData,2,segmentedRegression)
  #totalStats<-apply(newData,2,calculateCoefPow)
  #Calculate the average of the bootstrap statistics
  #averages<-apply(totalStats,1,mean)
  return(totalStats)
}

#finalResultsFre<-apply(data_fre[2:209],2,statBootstrap)
#finalResultsSpa<-apply(data_spa[2:209],2,statBootstrap)
#finalResultsGer<-apply(data_ger[2:209],2,statBootstrap)
#finalResultsIta<-apply(data_ita[2:209],2,statBootstrap)
#finalResultsEng<-apply(data_eng[2:209],2,statBootstrap)
#finalResultsRus<-apply(data_rus[2:209],2,statBootstrap)

finalResultsEng<-apply(data_eng[2:209],2,statBootstrap)
#finalResultsIta<-apply(data_ita[51:209],2,statBootstrap)
save(finalResultsEng,file="resultsEnglish")