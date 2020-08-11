#!/usr/bin/env Rscript
install.packages("igraph")
install.packages("MASS")
install.packages("stats")
install.packages("SOAR")
install.packages("segmented")
install.packages("psych")
install.packages("poweRlaw")
install.packages("psych")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("gtools")

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

#setwd("~/GoogleDrive/Project/Tables/")
setwd("")

data_ger = read.csv("all_words_per_year_ger.txt")
data_spa = read.csv("all_words_per_year_spa.txt")
data_fre = read.csv("all_words_per_year_fre.txt")
data_ita = read.csv("all_words_per_year_ita.txt")
data_eng = read.csv("all_words_per_year_eng.txt", quote = "", row.names = NULL, stringsAsFactors = FALSE)
data_rus = read.csv("all_words_per_year_rus.txt", quote = "", row.names = NULL, stringsAsFactors = FALSE)

data<-list(data_eng,data_ger,data_rus,data_spa,data_fre,data_ita)#,data_chi)
names(data)<-c("eng","ger","rus","spa","fre","ita")

isPercentage<-function(var,total){
  if(var<2){
    var<-floor(total*var)
  }
  return(var)
}

##Calculates Zipf's coefficient by power law fitting (Clauset's method) using all the data
#It returns the coefficient, the total number of tokens, two results concerning 
#goodness of fit tests and the value of x for which the power law starts fitting (The minimum frequency)

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
  total<-sum(as.numeric(occur_per_year))
  #occur_per_yearNorm<-occur_per_year/total
  #The sequence of ranks is created
  x <- seq_along(occur_per_year)
  xmin<-0
  #The coefficients of the linear regressions are calculated
  if(length(occur_per_year)>xmin)
  {  
    m <- power.law.fit(occur_per_year,implementation='plfit')
    #Normalized
    #m1<- power.law.fit(occur_per_yearNorm,implementation='plfit')
    value<-m$alpha
    #value1<-m1$alpha
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
  return(c(-value,total,test,test2,xmin))
  
}

##Calculates Zipf's coefficient by power law fitting (Clauset's method) 
##using only a fixed number of words.
#It returns the coefficient, the total number of tokens, two results concerning 
#goodness of fit tests and the value of x for which the power law starts fitting (The minimum frequency)
calculateCoefPowFix<-function(occur_per_year)
{
  #Possible NA values are filtered out
  occur_per_year<-occur_per_year[!is.na(occur_per_year)]
  #years with no occcurrence of a word are omitted, since the logarithm will not be defined.
  occur_per_year<-occur_per_year[occur_per_year!=0]
  #The words are ranked
  occur_per_year<-sort(occur_per_year, decreasing = TRUE)
  
  #Normalization
  total<-sum(as.numeric(occur_per_year))
  #occur_per_yearNorm<-occur_per_year/total
  #The sequence of ranks is created
  x <- seq_along(occur_per_year)
  xmin<-0
  #The coefficients of the linear regressions are calculated
  if(length(occur_per_year)>xmin)
  {  
    m <- power.law.fit(occur_per_year,implementation='plfit')
    #Normalized
    #m1<- power.law.fit(occur_per_yearNorm,implementation='plfit')
    value<-m$alpha
    #value1<-m1$alpha
    test<-m$KS.stat
    test2<-m$KS.p
    xmin<-m$xmin
  }
  else
  {
    value<-NA
    #value1<-NA
    test<-NA
    test2<-NA
    xmin<-NA
    
  }
  return(c(value,total,test,test2,xmin))
  
}

#results<-lapply(data,function(name){return(apply(name[2:209],2,calculateCoefPow,init=2,end=nrow(name)-1))})


##Calculates Zipf's coefficient for the kernel and unlimited lexicons using a segmented regression method
#It returns both coefficients, the proportion of words in the kernel lexicon and the value of the 
#transition point between the two different regimes.
segmentedRegression<-function(data)
{
  
  data<-data[!is.na(data)]
  data<-data[data!=0]
  totalData<-sum(as.numeric(data))
  maxDoc<-max(as.numeric(data))
  dataNorm<-data/totalData
  numWords<-length(data)
  data_hist<-try(hist(data,breaks = 1000000, plot=FALSE))
  if(class(data_hist) == "try-error") {return(c(NA,NA,NA,NA,NA))}
  x<-data_hist$mids[data_hist$counts!=0]#,decreasing=TRUE)
  y<-data_hist$counts[data_hist$counts!=0]/numWords#,decreasing=FALSE)
  
  x1<-x[x>1000]#x[3:(length(x)-2)]
  y<-y[x>1000]#y[3:(length(y)-2)]
  x<-x1
  #Where the power law stops holding
  aux<-power.law.fit(y,implementation="plfit")$xmin
  psi.init<-log(aux*totalData)
  #Plot empirical vs fitted distribution density functions
  x<-log(x)
  y<-log(y)
  
  
  lin.mod<-lm(y~x)
  #psi.init<-(log(aux)-lin.mod$coefficients[[1]])/lin.mod$coefficients[[2]]
  segmented.mod <- try(segmented(lin.mod, seg.Z = ~x, psi=15))
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
  return(c(floor(sum(as.numeric(b))),sum(as.numeric(c))/totalData,-coef1,-coef2,floor(exp(psi))))
}

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
  if(length(ranNums[is.na(ranNums)]>0)) print(i)
  #A data frame with the indexes of the words selected at random and the number of times they were selected
  indices<-as.data.frame(table(findInterval(ranNums,probs,all.inside=TRUE)))#sapply(ranNums,findPosition,cumVector=probs)))
  indices[[1]]<-as.numeric(as.character(indices[[1]]))
  newData<-rep(0,length(data))
  newData[indices[[1]]]<-indices[[2]]
  return(newData)
  
  
}

#Generates 100 bootstrapped samples and calculates Zipf's coefficients for them.
statBootstrap<-function(data)
{
  newData<-replicate(100,bootstrapData(data))
  totalStats<-apply(newData,2,segmentedRegression)
  #Calculate the average of the bootstrap statistics
  #averages<-apply(totalStats,1,mean)
  #print("aaaaaa")
  return(totalStats)
}

#Applies segmented regression over all the data
segReg<-lapply(data,function(name){return(apply(name[2:209],2,segmentedRegression))})

#segEng<-apply(data_eng[2:209],2,segmentedRegression)
#segGer<-apply(data_ger[2:209],2,segmentedRegression)
#segRus<-apply(data_rus[2:209],2,segmentedRegression)
#segSpa<-apply(data_spa[2:209],2,segmentedRegression)
#segFre<-apply(data_fre[2:209],2,segmentedRegression)
#segIta<-apply(data_ita[2:209],2,segmentedRegression)

#Separates the results according to languages
fixKer<-lapply(segReg,function(name){return(name[1,])})
names(fixKer)<-c("eng","ger","rus","spa","fre","ita")
fixKerRel<-lapply(segReg,function(name){return(name[2,])})
names(fixKerRel)<-c("eng","ger","rus","spa","fre","ita")
alpha1<-lapply(segReg,function(name){return(name[3,])})
names(alpha1)<-c("eng","ger","rus","spa","fre","ita")
alpha2<-lapply(segReg,function(name){return(name[4,])})
names(alpha2)<-c("eng","ger","rus","spa","fre","ita")
kink<-lapply(segReg,function(name){return(name[5,])})
names(kink)<-c("eng","ger","rus","spa","fre","ita")

##Cleaning the final data to produce meaningful plots
#Some exceptions are eliminated (the points won't appeared in the plot)
alpha1$spa[alpha1$spa>1.2 | alpha1$spa<0.8]<-NA
alpha1$ita[alpha1$ita>1.15 | alpha1$ita<0.8]<-NA
alpha1$ger[alpha1$ger>1.2 | alpha1$ger<0.95]<-NA
alpha1$rus[alpha1$rus>1.2 | alpha1$rus<0.8]<-NA
alpha2$spa[alpha2$spa>4]<-NA
alpha2$rus[alpha2$rus>4]<-NA
alpha2$ita[alpha2$ita>4]<-NA
alpha2$ger[alpha2$ger>4]<-NA
alpha2$eng[alpha2$eng>4]<-NA
kink$spa[kink$spa>10000]<-NA
kink$eng[kink$eng>10000]<-NA

#Plots the results for all languages in one single figure.
dataFrameEng<-data.frame(year,fixKerRel$eng,rep("English",208),row.names=NULL)
names(dataFrameEng)<-c("year","coef","lang")
dataFrameGer<-data.frame(year,fixKerRel$ger,rep("German",208),row.names=NULL)
names(dataFrameGer)<-c("year","coef","lang")
dataFrameRus<-data.frame(year,fixKerRel$rus,rep("Russian",208),row.names=NULL)
names(dataFrameRus)<-c("year","coef","lang")
dataFrameSpa<-data.frame(year,fixKerRel$spa,rep("Spanish",208),row.names=NULL)
names(dataFrameSpa)<-c("year","coef","lang")
dataFrameFre<-data.frame(year,fixKerRel$fre,rep("French",208),row.names=NULL)
names(dataFrameFre)<-c("year","coef","lang")
dataFrameIta<-data.frame(year,fixKerRel$ita,rep("Italian",208),row.names=NULL)
names(dataFrameIta)<-c("year","coef","lang")


df<-rbind(dataFrameEng,dataFrameGer,dataFrameRus,dataFrameSpa,dataFrameFre,dataFrameIta)

g<-ggplot(df, aes(year, coef)) + geom_line(color="coral3")+
  facet_wrap( ~ lang, ncol = 2,scales="free_y") +
  #theme_bw() + opts(strip.background=theme_blank())
  labs(x="Year", y=expression(alpha[2]))+#ylim(0.8,1.2)+
  #geom_abline(intercept = lm(varname~year)$coefficients[[1]], lm(varname~year)$coefficients[[2]],color="blue")+
  theme_light()+
  theme(strip.text = element_text(size=20),plot.title = element_text(size=18, face="bold", vjust=2),axis.title.x = element_text(size=20,color="black", vjust=-0.5),
        axis.title.y = element_text(size=20,color="black" , vjust=0.5) )

g

#######
#######

##Generates bootstarps for all languages
finalResultsFre<-apply(data_fre[2:209],2,statBootstrap)
finalResultsSpa<-apply(data_spa[2:209],2,statBootstrap)
finalResultsGer<-apply(data_ger[2:209],2,statBootstrap)
finalResultsIta<-apply(data_ita[2:209],2,statBootstrap)
finalResultsEng<-apply(data_eng[2:209],2,statBootstrap)
finalResultsRus<-apply(data_rus[2:209],2,statBootstrap)


#Calculates Zipf's coefficient of the 5000 highest frequency words. It is done via a 
#linear regresion because we want to check that the evolution observed in the kernel lexicon
#is not due to size, and using Clauset's method the number of words taken into account would 
#change for each bootstrapped sample. Less words can be considered (e.g. 2000) to discard any 
#possible influence of the tail.
calculateCoefFix<-function(data)
{
  data<-data[data!=0]
  data<-data[!is.na(data)]
  totalData<-sum(as.numeric(data))
  maxDoc<-max(data)
 
  data<-data/totalData
  numWords<-length(data)
  data_hist<-try(hist(log(data),breaks = 30, plot=FALSE))
  if(class(data_hist) == "try-error") {return(NA)}
  x<-sort(data_hist$counts[data_hist$counts!=0],decreasing=FALSE)
  y<-sort(exp(data_hist$mids[data_hist$counts!=0]),decreasing=TRUE)
  x<-x[1:(length(x)-2)]
  y<-y[1:(length(y)-2)]
  
  #Plot empirical vs fitted distribution density functions
  x<-log(x)
  y<-log(y)
  
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

#Selects the words with ordered according to frequency in the range [5000,50000] .
fix2<-function(data)
{
  data<-sort(data,decreasing=TRUE)
  return(data[5000:50000])
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

##Plots the results for one language (french in this case)
finalResultsFreFix<-apply(data_fre[2:209],2,statBootstrapFix)
#finalResultsSpaFix<-apply(data_spa[2:209],2,statBootstrapFix)
#finalResultsEngFix<-apply(data_eng[2:209],2,statBootstrapFix)

#boxplot(finalResultsSpaFix,outline=FALSE,xlab="year",ylim=c(0.95,1.1),ylab=expression(alpha[1]),main=expression(paste("Evolution of ", alpha[1], " for Spanish")))
finFre<-melt(as.data.frame(finalResultsFreFix))
finFre[,1]<-as.numeric(finFre[,1])+1799
ggplot(finFre, aes(x=variable, y=value,group=variable)) + geom_boxplot(color="firebrick")+
  ggtitle(expression(paste('Evolution of ', alpha,' for the fixed kernel lexicon in French')))+
  labs(x="Year", y=expression(alpha))+
  theme(plot.title = element_text(size=20, face="bold", vjust=2),axis.title.x = element_text(size=15, color="black", vjust=-0.5),axis.title.y = element_text(size=15,color="black" , vjust=0.5)) 

separatingResults<-function(data,n)
{
  return(data[seq(1,600,1)%%6==n])
}

#####
##Experimental plots (Just to check results, not suitable for the report)
#####

load("Spa.RDa")
load("Eng.RData")
load("Rus.RDa")
load("Ita.RDa")
load("Ger.RDa")
load("French.RData")


boxplot(relativeKerGer,outline=FALSE,xlab="year",ylab="Percentage of vocabulary in the kernel lexicon",main="Evolution of the relative size of the kernel lexicon for German",names=as.character(seq(1800,2007,1)))
boxplot(alpha1Ger,outline=FALSE,xlab="year",ylim=c(0.9,1.1),ylab="Alpha1",main="Evolution of alpha1 for German")
boxplot(alpha2Ger,outline=FALSE,ylim=c(1,3),xlab="year",ylab="Alpha2",main="Evolution of alpha2 for German")
boxplot(xminGer,outline=FALSE,xlab="year",ylab="xmin",main="Evolution of the position until which the power law holds for German")
boxplot(kinkGer,outline=FALSE,xlab="year",ylim=c(3,8),ylab="Kink",main="Evolution of the kink for German")

####################################
#Plotting the final results
####################################

#Cleans the data removing NA's
finalResultsRus<-as.data.frame(finalResultsRus)
for(i in seq(0,17,1))
{name=paste("X",1800+i)
 name<-gsub(" ","",name)
  finalResultsRus[[name]]<-rep(NA,600)}
finalResultsIta<-as.data.frame(finalResultsIta)
for(i in seq(0,1,1))
{name=paste("X",1800+i)
 name<-gsub(" ","",name)
 finalResultsIta[[name]]<-rep(NA,600)}


##Separates the different variables returned by the functions calculating Zipf's coefficients
relativeKerEng<-as.data.frame(apply(finalResultsEng,2,separatingResults,n=1))
alpha1Eng<-as.data.frame(apply(finalResultsEng,2,separatingResults,n=2))
alpha2Eng<-as.data.frame(apply(finalResultsEng,2,separatingResults,n=3))
xminEng<-as.data.frame(apply(finalResultsEng,2,separatingResults,n=4))
kinkEng<-as.data.frame(apply(finalResultsEng,2,separatingResults,n=5))
psi.initEng<-as.data.frame(apply(finalResultsEng,2,separatingResults,n=0))

relativeKerGer<-as.data.frame(apply(finalResultsGer,2,separatingResults,n=1))
alpha1Ger<-as.data.frame(apply(finalResultsGer,2,separatingResults,n=2))
alpha2Ger<-as.data.frame(apply(finalResultsGer,2,separatingResults,n=3))
xminGer<-as.data.frame(apply(finalResultsGer,2,separatingResults,n=4))
kinkGer<-as.data.frame(apply(finalResultsGer,2,separatingResults,n=5))
kinkGer<-as.data.frame(apply(finalResultsGer,2,separatingResults,n=5))

#finalResultsRus<-lapply(finalResultsGer,function(name){return(as.vector(unlist(name,use.names=FALSE)))})
relativeKerRus<-as.data.frame(apply(finalResultsRus,2,separatingResults,n=1))
alpha1Rus<-as.data.frame(apply(finalResultsRus,2,separatingResults,n=2))
alpha2Rus<-as.data.frame(apply(finalResultsRus,2,separatingResults,n=3))
xminRus<-as.data.frame(apply(finalResultsRus,2,separatingResults,n=4))
kinkRus<-as.data.frame(apply(finalResultsRus,2,separatingResults,n=5))
psi.initRus<-as.data.frame(apply(finalResultsRus,2,separatingResults,n=0))

relativeKerSpa<-as.data.frame(apply(finalResultsSpa,2,separatingResults,n=1))
alpha1Spa<-as.data.frame(apply(finalResultsSpa,2,separatingResults,n=2))
alpha2Spa<-as.data.frame(apply(finalResultsSpa,2,separatingResults,n=3))
xminSpa<-as.data.frame(apply(finalResultsSpa,2,separatingResults,n=4))
kinkSpa<-as.data.frame(apply(finalResultsSpa,2,separatingResults,n=5))
psi.initSpa<-as.data.frame(apply(finalResultsSpa,2,separatingResults,n=0))

relativeKerFre<-as.data.frame(apply(finalResultsFre,2,separatingResults,n=1))
alpha1Fre<-as.data.frame(apply(finalResultsFre,2,separatingResults,n=2))
alpha2Fre<-as.data.frame(apply(finalResultsFre,2,separatingResults,n=3))
xminFre<-as.data.frame(apply(finalResultsFre,2,separatingResults,n=4))
kinkFre<-as.data.frame(apply(finalResultsFre,2,separatingResults,n=5))
psi.initFre<-as.data.frame(apply(finalResultsFre,2,separatingResults,n=0))

relativeKerIta<-as.data.frame(apply(finalResultsIta,2,separatingResults,n=1))
alpha1Ita<-as.data.frame(apply(finalResultsIta,2,separatingResults,n=2))
alpha2Ita<-as.data.frame(apply(finalResultsIta,2,separatingResults,n=3))
xminIta<-as.data.frame(apply(finalResultsIta,2,separatingResults,n=4))
kinkIta<-as.data.frame(apply(finalResultsIta,2,separatingResults,n=5))
psi.initIta<-as.data.frame(apply(finalResultsIta,2,separatingResults,n=0))


##Plotting the bootstrapped kernel for one language (example with english)

a2Eng<-melt(as.data.frame(alpha2Eng))
a2Eng[,1]<-as.numeric(a2Eng[,1])+1799
ggplot(a2Eng, aes(x=variable, y=value,group=variable)) + geom_boxplot(color="firebrick")+
  ggtitle(expression(paste('Evolution of ', alpha[2],' for English')))+
  labs(x="Year", y=expression(alpha[2]))+ylim(1.5,2.2)
  theme(plot.title = element_text(size=20, face="bold", vjust=2),axis.title.x = element_text(size=20, color="black", vjust=-0.5),axis.title.y = element_text(size=20,color="black" , vjust=0.5)) 
####
#Plotting bootstraps
####

###
#Alpha2
###
A2E<-melt(as.data.frame(relativeKerEng))
A2G<-melt(as.data.frame(relativeKerGer))
A2R<-melt(as.data.frame(relativeKerRus))
A2S<-melt(as.data.frame(relativeKerSpa))
A2F<-melt(as.data.frame(relativeKerFre))
A2I<-melt(as.data.frame(relativeKerIta))


A2E[,1]<-as.numeric(A2E[,1])+1799
A2G[,1]<-as.numeric(A2G[,1])+1799
A2R[,1]<-as.numeric(A2R[,1])+1799
A2S[,1]<-as.numeric(A2S[,1])+1799
A2F[,1]<-as.numeric(A2F[,1])+1799
A2I[,1]<-as.numeric(A2I[,1])+1799

A2E$value[A2E$value>4]<-NA
A2G$value[A2G$value>4]<-NA
A2R$value[A2R$value>4]<-NA
A2S$value[A2S$value>4]<-NA
A2F$value[A2F$value>4]<-NA
A2I$value[A2I$value>4]<-NA

ggplot(A2E, aes(x=variable, y=value,group=variable)) + geom_boxplot(color="firebrick")+
  ggtitle(expression(paste('Evolution of ', alpha[1],' for English')))+
  labs(x="Year", y=expression(alpha[2]))+
  theme(plot.title = element_text(size=20, face="bold", vjust=2),axis.title.x = element_text(size=20, color="black", vjust=-0.5),axis.title.y = element_text(size=20,color="black" , vjust=0.5)) 

A2F<-melt(as.data.frame(BootsFreFix))
A2F[,1]<-as.numeric(A2F[,1])+1799
ggplot(BFF, aes(x=variable, y=value,group=variable)) + geom_boxplot(color="firebrick")+
  ggtitle(expression(paste('Evolution of ', alpha[1],' for French')))+
  labs(x="Year", y=expression(alpha[2]))+
  theme(plot.title = element_text(size=20, face="bold", vjust=2),axis.title.x = element_text(size=20, color="black", vjust=-0.5),axis.title.y = element_text(size=20,color="black" , vjust=0.5)) 


dataFrameEng<-data.frame(A2E,rep("English",2080),row.names=NULL)
names(dataFrameEng)<-c("year","coef","lang")
dataFrameGer<-data.frame(A2G,rep("German",2080),row.names=NULL)
names(dataFrameGer)<-c("year","coef","lang")
dataFrameRus<-data.frame(A2R,rep("Russian",2080),row.names=NULL)
names(dataFrameRus)<-c("year","coef","lang")
dataFrameSpa<-data.frame(A2S,rep("Spanish",2080),row.names=NULL)
names(dataFrameSpa)<-c("year","coef","lang")
dataFrameFre<-data.frame(A2F,rep("French",2080),row.names=NULL)
names(dataFrameFre)<-c("year","coef","lang")
dataFrameIta<-data.frame(A2I,rep("Italian",2080),row.names=NULL)
names(dataFrameIta)<-c("year","coef","lang")


df<-rbind(dataFrameEng,dataFrameGer,dataFrameRus,dataFrameSpa,dataFrameFre,dataFrameIta)

g<-ggplot(df, aes(year, coef, group=year)) + geom_boxplot(color="indianred3",outlier.colour = "black", outlier.size = 0.8)+
  facet_wrap( ~ lang, ncol = 2,scales="free_y") +
  #theme_bw() + opts(strip.background=theme_blank())
  labs(x="Year", y=expression(alpha[1]))+#ylim(1.95,y)+
  #geom_abline(intercept = lm(varname~year)$coefficients[[1]], lm(varname~year)$coefficients[[2]],color="blue")+
  theme_light()+
  theme(strip.text = element_text(size=20),plot.title = element_text(size=18, face="bold", vjust=2),axis.title.x = element_text(size=20,color="black", vjust=-0.5),
        axis.title.y = element_text(size=20,color="black" , vjust=0.5) )
g

####
#Plotting results for Power law fitting without bootstrap
#Changing the variables different quantities such as the coefficients or the location of the xmin
#can be plotted (e.g. instead of creating the data frames below with coefs$eng, use xmin$eng,
#and make the respective changes in the plotting function)
####
results<-lapply(data,function(name){return(apply(name[2:209],2,calculateCoefPow,init=2,end=nrow(name)-1))})
names(results)<-c("ger","spa","fre","ita","eng","rus")
total<-lapply(results,function(name){return(name[3,])})
names(total)<-c("ger","spa","fre","ita","eng","rus")
coefs<-lapply(results,function(name){return(-name[1,])})
names(coefs)<-c("ger","spa","fre","ita","eng","rus")
coefsNorm<-lapply(results,function(name){return(name[2,])})
names(coefsNorm)<-c("ger","spa","fre","ita","eng","rus")
xmin<-lapply(results,function(name){return(name[5,])})
names(xmin)<-c("ger","spa","fre","ita","eng","rus")


year<-seq(1800,2007,1)
plot(year,-coefs$eng,col="black",type="l",xlab="Year",ylab="Zipf's coefficient",main="Evolution of Zipf's coefficient for kernel lexicon")
plot(year,coefsFix$eng,col="black",type="l",xlab="Year",ylab="Zipf's coefficient",main="Evolution of Zipf's coefficient for kernel lexicon")


coefs$rus[1:25]<-NA
dataFrameEng<-data.frame(year,coefs$eng,rep("English",208),row.names=NULL)
names(dataFrameEng)<-c("year","coef","lang")
dataFrameGer<-data.frame(year,coefs$ger,rep("German",208),row.names=NULL)
names(dataFrameGer)<-c("year","coef","lang")
dataFrameRus<-data.frame(year,coefs$rus,rep("Russian",208),row.names=NULL)
names(dataFrameRus)<-c("year","coef","lang")
dataFrameSpa<-data.frame(year,coefs$spa,rep("Spanish",208),row.names=NULL)
names(dataFrameSpa)<-c("year","coef","lang")
dataFrameFre<-data.frame(year,coefs$fre,rep("French",208),row.names=NULL)
names(dataFrameFre)<-c("year","coef","lang")
dataFrameIta<-data.frame(year,coefs$ita,rep("Italian",208),row.names=NULL)
names(dataFrameIta)<-c("year","coef","lang")


df<-rbind(dataFrameEng,dataFrameGer,dataFrameRus,dataFrameSpa,dataFrameFre,dataFrameIta)

g<-ggplot(df, aes(year, coef)) + geom_line(color="skyblue3")+
  facet_wrap( ~ lang, ncol = 2,scales="free_y") +
  #theme_bw() + opts(strip.background=theme_blank())
  labs(x="Year", y=expression(alpha[2])))+#ylim(1.95,y)+
  #geom_abline(intercept = lm(varname~year)$coefficients[[1]], lm(varname~year)$coefficients[[2]],color="blue")+
  theme_light()+
  theme(strip.text = element_text(size=20),plot.title = element_text(size=18, face="bold", vjust=2),axis.title.x = element_text(size=20,color="black", vjust=-0.5),
        axis.title.y = element_text(size=20,color="black" , vjust=0.5) )

g

####
##Plotting fix coefficient and coefficient of the unlimited lexicon under bootstrap
####

coefs<-lapply(data,function(name){return(apply(name[2:209],2,fix2))})
names(coefs)<-c("eng","ger","rus","spa","fre","ita")
resultsFix<-lapply(coefs,function(name){return(apply(name,2,calculateCoefFix))})
resultsFixEng<-apply(dataFixKer$eng[2:209],2,calculateCoefPowFix)
names(resultsFix)<-c("eng","ger","rus","spa","fre","ita")
totalFix<-lapply(resultsFix,function(name){return(name[2,])})
names(totalFix)<-c("eng","ger","rus","spa","fre","ita")
coefsFix<-lapply(resultsFix,function(name){return(name[1,])})
names(coefsFix)<-c("eng","ger","rus","spa","fre","ita")
testFix<-lapply(resultsFix,function(name){return(name[3,])})
names(testFix)<-c("eng","ger","rus","spa","fre","ita")
xminFix<-lapply(resultsFix,function(name){return(name[5,])})
names(xminFix)<-c("eng","ger","rus","spa","fre","ita")

coefs<-resultsFix

coefsFix$rus[1:25]<-NA
dataFrameEng<-data.frame(year,coefs$eng,rep("English",208),row.names=NULL)
names(dataFrameEng)<-c("year","coef","lang")
dataFrameGer<-data.frame(year,coefs$ger,rep("German",208),row.names=NULL)
names(dataFrameGer)<-c("year","coef","lang")
dataFrameRus<-data.frame(year,coefs$rus,rep("Russian",208),row.names=NULL)
names(dataFrameRus)<-c("year","coef","lang")
dataFrameSpa<-data.frame(year,coefs$spa,rep("Spanish",208),row.names=NULL)
names(dataFrameSpa)<-c("year","coef","lang")
dataFrameFre<-data.frame(year,coefs$fre,rep("French",208),row.names=NULL)
names(dataFrameFre)<-c("year","coef","lang")
dataFrameIta<-data.frame(year,coefs$ita,rep("Italian",208),row.names=NULL)
names(dataFrameIta)<-c("year","coef","lang")


df<-rbind(dataFrameEng,dataFrameGer,dataFrameRus,dataFrameSpa,dataFrameFre,dataFrameIta)

g<-ggplot(df, aes(year, coef)) + geom_line(color="indianred3")+
  facet_wrap( ~ lang, ncol = 2,scales="free_y") +
  #theme_bw() + opts(strip.background=theme_blank())
  labs(x="Year", y=expression(alpha[2]))+#ylim(1.95,y)+
  #geom_abline(intercept = lm(varname~year)$coefficients[[1]], lm(varname~year)$coefficients[[2]],color="blue")+
  theme_light()+
  theme(strip.text = element_text(size=20),plot.title = element_text(size=18, face="bold", vjust=2),axis.title.x = element_text(size=20,color="black", vjust=-0.5),
        axis.title.y = element_text(size=20,color="black" , vjust=0.5) )

g


plotDataFrame<-data.frame(year,coefsFix$ita)
ggplot(plotDataFrame, aes(year, coefsFix.ita))+geom_line(color="firebrick")+
  ggtitle("Evolution of Zipf's coefficient for English")+
  labs(x="Year", y="Zipf's coefficient")+#ylim(1.8,2.3)+
  theme(plot.title = element_text(size=20, face="bold", vjust=2),axis.title.x = element_text(size=20,color="black", vjust=-0.5),
        axis.title.y = element_text(size=20,color="black" , vjust=0.5) )

################
###Plotting results for fixed kernel (5000 words)
################

#Macro to generate the fixed kernel plots for all languages 
#Another way of plotting, taken from http://www.r-bloggers.com/producing-grids-of-plots-in-r-with-ggplot2-a-journey-of-discovery/
#Didn't completely work
cg <- defmacro(dataFrame,varname, vartext,y, expr={ggplot(dataFrame, aes(year, varname)) + geom_line(color="firebrick")+
            ggtitle(vartext)+
              labs(x="Year", y="Zipf's coefficient")+ylim(1.95,y)+
              #geom_abline(intercept = lm(varname~year)$coefficients[[1]], lm(varname~year)$coefficients[[2]],color="blue")+
              theme(plot.title = element_text(size=18, face="bold", vjust=2),axis.title.x = element_text(size=15,color="black", vjust=-0.5),
                    axis.title.y = element_text(size=15,color="black" , vjust=0.5) )})

#Classical method
# Create all of the graphs we want, storing them in variables

eng = cg(data.frame(year,coefsFix$eng), coefsFix$eng,"English",2.2)
ger = cg(data.frame(year,coefsFix$ger), coefsFix$ger,"German",2.11)
rus = cg(data.frame(year,coefsFix$rus), coefsFix$rus,"Russian",2.2)
spa = cg(data.frame(year,coefsFix$spa), coefsFix$spa,"Spanish",2.2)
fra = cg(data.frame(year,coefsFix$fre), coefsFix$fre,"French",2.1)
ita = cg(data.frame(year,coefsFix$ita), coefsFix$ita,"Italian",2.1)

arrange_ggplot2(eng,ger,rus,spa,fra,ita, ncol=2)


coefsFix$rus[1:25]<-NA
dataFrameEng<-data.frame(year,coefsFix$eng,rep("English",208),row.names=NULL)
names(dataFrameEng)<-c("year","coef","lang")
dataFrameGer<-data.frame(year,coefsFix$ger,rep("German",208),row.names=NULL)
names(dataFrameGer)<-c("year","coef","lang")
dataFrameRus<-data.frame(year,coefsFix$rus,rep("Russian",208),row.names=NULL)
names(dataFrameRus)<-c("year","coef","lang")
dataFrameSpa<-data.frame(year,coefsFix$spa,rep("Spanish",208),row.names=NULL)
names(dataFrameSpa)<-c("year","coef","lang")
dataFrameFre<-data.frame(year,coefsFix$fre,rep("French",208),row.names=NULL)
names(dataFrameFre)<-c("year","coef","lang")
dataFrameIta<-data.frame(year,coefsFix$ita,rep("Italian",208),row.names=NULL)
names(dataFrameIta)<-c("year","coef","lang")


df<-rbind(dataFrameEng,dataFrameGer,dataFrameRus,dataFrameSpa,dataFrameFre,dataFrameIta)

g<-ggplot(df, aes(year, coef)) + geom_line(color="firebrick")+
  facet_wrap( ~ lang, ncol = 2,scales="free_y") +
  #theme_bw() + opts(strip.background=theme_blank())
  labs(x="Year", y="Zipf's coefficient")+#ylim(1.95,y)+
  #geom_abline(intercept = lm(varname~year)$coefficients[[1]], lm(varname~year)$coefficients[[2]],color="blue")+
  theme_light()+
  theme(strip.text = element_text(size=20),plot.title = element_text(size=18, face="bold", vjust=2),axis.title.x = element_text(size=20,color="black", vjust=-0.5),
        axis.title.y = element_text(size=20,color="black" , vjust=0.5) )

g

####
##Plotting original coefficient (not fixed kernel) without bootstrap
####
coefsFix$rus[1:25]<-NA
dataFrameEng<-data.frame(year,coefs$eng,rep("English",208),row.names=NULL)
names(dataFrameEng)<-c("year","coef","lang")
dataFrameGer<-data.frame(year,coefs$ger,rep("German",208),row.names=NULL)
names(dataFrameGer)<-c("year","coef","lang")
dataFrameRus<-data.frame(year,coefs$rus,rep("Russian",208),row.names=NULL)
names(dataFrameRus)<-c("year","coef","lang")
dataFrameSpa<-data.frame(year,coefs$spa,rep("Spanish",208),row.names=NULL)
names(dataFrameSpa)<-c("year","coef","lang")
dataFrameFre<-data.frame(year,coefs$fre,rep("French",208),row.names=NULL)
names(dataFrameFre)<-c("year","coef","lang")
dataFrameIta<-data.frame(year,coefs$ita,rep("Italian",208),row.names=NULL)
names(dataFrameIta)<-c("year","coef","lang")


df<-rbind(dataFrameEng,dataFrameGer,dataFrameRus,dataFrameSpa,dataFrameFre,dataFrameIta)

g<-ggplot(df, aes(year, coef)) + geom_line(color="skyblue3")+
  facet_wrap( ~ lang, ncol = 2,scales="free_y") +
  #theme_bw() + opts(strip.background=theme_blank())
  labs(x="Year", y="Zipf's coefficient")+#ylim(1.95,y)+
  #geom_abline(intercept = lm(varname~year)$coefficients[[1]], lm(varname~year)$coefficients[[2]],color="blue")+
  theme_light()+
  theme(strip.text = element_text(size=20),plot.title = element_text(size=18, face="bold", vjust=2),axis.title.x = element_text(size=20,color="black", vjust=-0.5),
        axis.title.y = element_text(size=20,color="black" , vjust=0.5) )

g

####
##Plotting relative size xmin
####

findKerLim<-function(data)
{
  xmin<-data[length(data)]
  data<-data[1:(length(data)-1)]
  data<-data[data!=0]
  absker<-length(data[data>xmin])
  return(c(absker/length(data),absker))
}

dataEng<-rbind(data_eng[,2:209],xmin$eng)
kerEng<-apply(dataEng,2,findKerLim)
relEng<-kerEng[1,]
absEng<-kerEng[2,]
dataGer<-rbind(data_ger[,2:209],xmin$ger)
kerGer<-apply(dataGer,2,findKerLim)
relGer<-kerGer[1,]
absGer<-kerGer[2,]
dataRus<-rbind(data_rus[,2:209],xmin$rus)
kerRus<-apply(dataSpa,2,findKerLim)
relRus<-kerRus[1,]
absRus<-kerRus[2,]
dataSpa<-rbind(data_spa[,2:209],xmin$spa)
kerSpa<-apply(dataSpa,2,findKerLim)
relSpa<-kerSpa[1,]
absSpa<-kerSpa[2,]
dataFre<-rbind(data_fre[,2:209],xmin$fre)
kerFre<-apply(dataFre,2,findKerLim)
relFre<-kerFre[1,]
absFre<-kerFre[2,]
dataIta<-rbind(data_ita[,2:209],xmin$ita)
kerIta<-apply(dataIta,2,findKerLim)
relIta<-kerIta[1,]
absIta<-kerIta[2,]

xmin$rus[1:25]<-NA
dataFrameEng<-data.frame(year,relEng,rep("English",208),row.names=NULL)
names(dataFrameEng)<-c("year","coef","lang")
dataFrameGer<-data.frame(year,relGer,rep("German",208),row.names=NULL)
names(dataFrameGer)<-c("year","coef","lang")
dataFrameRus<-data.frame(year,relRus,rep("Russian",208),row.names=NULL)
names(dataFrameRus)<-c("year","coef","lang")
dataFrameSpa<-data.frame(year,relSpa,rep("Spanish",208),row.names=NULL)
names(dataFrameSpa)<-c("year","coef","lang")
dataFrameFre<-data.frame(year,relFre,rep("French",208),row.names=NULL)
names(dataFrameFre)<-c("year","coef","lang")
dataFrameIta<-data.frame(year,relIta,rep("Italian",208),row.names=NULL)
names(dataFrameIta)<-c("year","coef","lang")


dataFrameEng<-data.frame(year,absEng,rep("English",208),row.names=NULL)
names(dataFrameEng)<-c("year","coef","lang")
dataFrameGer<-data.frame(year,absGer,rep("German",208),row.names=NULL)
names(dataFrameGer)<-c("year","coef","lang")
dataFrameRus<-data.frame(year,absRus,rep("Russian",208),row.names=NULL)
names(dataFrameRus)<-c("year","coef","lang")
dataFrameSpa<-data.frame(year,absSpa,rep("Spanish",208),row.names=NULL)
names(dataFrameSpa)<-c("year","coef","lang")
dataFrameFre<-data.frame(year,absFre,rep("French",208),row.names=NULL)
names(dataFrameFre)<-c("year","coef","lang")
dataFrameIta<-data.frame(year,absIta,rep("Italian",208),row.names=NULL)
names(dataFrameIta)<-c("year","coef","lang")

dataFrameGer$coef[dataFrameGer$coef>5000]<-NA
df<-rbind(dataFrameEng,dataFrameGer,dataFrameRus,dataFrameSpa,dataFrameFre,dataFrameIta)

g<-ggplot(df, aes(year, coef)) + geom_line(color="skyblue3")+
  facet_wrap( ~ lang, ncol = 2,scales="free_y") +
  #theme_bw() + opts(strip.background=theme_blank())
  labs(x="Year", y="xmin")+#ylim(1.95,y)+
  #geom_abline(intercept = lm(varname~year)$coefficients[[1]], lm(varname~year)$coefficients[[2]],color="blue")+
  theme_light()+
  theme(strip.text = element_text(size=20),plot.title = element_text(size=18, face="bold", vjust=2),axis.title.x = element_text(size=20,color="black", vjust=-0.5),
        axis.title.y = element_text(size=20,color="black" , vjust=0.5) )

g

###########
#Plotting bootstrapped results for the fixed kernel
###########

load("BootsEngFix2.rda")
load("BootsGerFix2.rda")
load("BootsRusFix2.rda")
load("BootsSpaFix2.rda")
load("BootsFreFix2.rda")
load("BootsItaFix2.rda")

#BootsEngFix<-apply(data_eng[2:209],2,statBootstrapFix)
#BootsGerFix<-apply(data_ger[2:209],2,statBootstrapFix)
#BootsRusFix<-apply(data_rus[2:209],2,statBootstrapFix)
#BootsSpaFix<-apply(data_spa[2:209],2,statBootstrapFix)
#BootsFreFix<-apply(data_fre[2:209],2,statBootstrapFix)
#BootsItaFix<-apply(data_ita[2:209],2,statBootstrapFix)

BootsRusFix<-as.data.frame(BootsRusFix)
for(i in seq(0,17,1))
{name=paste("X",1800+i)
 name<-gsub(" ","",name)
 BootsRusFix[[name]]<-rep(NA,100)}
BootsRusFix<-BootsRusFix[,order(names(BootsRusFix))]
BootsRusFix<-as.matrix(BootsRusFix)

BootsItaFix$X1802<-NULL
BootsItaFix$X1810<-NULL
BootsItaFix$X1812<-NULL
BootsItaFix$X1814<-NULL
BootsItaFix<-as.data.frame(BootsItaFix)
for(i in seq(0,14,1))
{name=paste("X",1800+i)
 name<-gsub(" ","",name)
 if(is.null(BootsItaFix[[name]])) BootsItaFix[[name]]<-rep(NA,100)}
BootsItaFix<-BootsItaFix[,order(names(BootsItaFix))]
BootsItaFix<-as.matrix(BootsItaFix)

BEF<-melt(as.data.frame(BootsEngFix))
BGF<-melt(as.data.frame(BootsGerFix))
BRF<-melt(as.data.frame(BootsRusFix))
BSF<-melt(as.data.frame(BootsSpaFix))
BFF<-melt(as.data.frame(BootsFreFix))
BIF<-melt(as.data.frame(BootsItaFix))

BEF[,1]<-as.numeric(BEF[,1])+1799
BGF[,1]<-as.numeric(BGF[,1])+1799
BRF[,1]<-as.numeric(BRF[,1])+1799
BSF[,1]<-as.numeric(BSF[,1])+1799
BFF[,1]<-as.numeric(BFF[,1])+1799
BIF[,1]<-as.numeric(BIF[,1])+1799

BEF$value[BEF$value>1.15]<-NA
BGF$value[BGF$value>1.15]<-NA
BRF$value[BRF$value>1.1]<-NA
BSF$value[BSF$value>1.1]<-NA
BFF$value[BFF$value>1.15]<-NA
BIF$value[BIF$value>1.1]<-NA


##Examples for english and french
ggplot(BEF, aes(x=variable, y=value,group=variable)) + geom_boxplot(color="firebrick")+
  ggtitle(expression(paste('Evolution of ', alpha[1],' for English')))+
  labs(x="Year", y=expression(alpha[1]))+
theme(plot.title = element_text(size=20, face="bold", vjust=2),axis.title.x = element_text(size=20, color="black", vjust=-0.5),axis.title.y = element_text(size=20,color="black" , vjust=0.5)) 



ggplot(BFF, aes(x=variable, y=value,group=variable)) + geom_boxplot(color="firebrick")+
  ggtitle(expression(paste('Evolution of ', alpha[1],' for French')))+
  labs(x="Year", y=expression(alpha[1]))+
  theme(plot.title = element_text(size=20, face="bold", vjust=2),axis.title.x = element_text(size=20, color="black", vjust=-0.5),axis.title.y = element_text(size=20,color="black" , vjust=0.5)) 


dataFrameEng<-data.frame(BEF,rep("English",2080),row.names=NULL)
names(dataFrameEng)<-c("year","coef","lang")
dataFrameGer<-data.frame(BGF,rep("German",2080),row.names=NULL)
names(dataFrameGer)<-c("year","coef","lang")
dataFrameRus<-data.frame(BRF,rep("Russian",2080),row.names=NULL)
names(dataFrameRus)<-c("year","coef","lang")
dataFrameSpa<-data.frame(BSF,rep("Spanish",2080),row.names=NULL)
names(dataFrameSpa)<-c("year","coef","lang")
dataFrameFre<-data.frame(BFF,rep("French",2080),row.names=NULL)
names(dataFrameFre)<-c("year","coef","lang")
dataFrameIta<-data.frame(BIF,rep("Italian",2080),row.names=NULL)
names(dataFrameIta)<-c("year","coef","lang")


df<-rbind(dataFrameEng,dataFrameGer,dataFrameRus,dataFrameSpa,dataFrameFre,dataFrameIta)

g<-ggplot(df, aes(year, coef, group=year)) + geom_boxplot(color="indianred3",outlier.colour = "black", outlier.size = 0.8)+
  facet_wrap( ~ lang, ncol = 2,scales="free_y") +
  #theme_bw() + opts(strip.background=theme_blank())
  labs(x="Year", y=expression(alpha[1]))+#ylim(1.95,y)+
  #geom_abline(intercept = lm(varname~year)$coefficients[[1]], lm(varname~year)$coefficients[[2]],color="blue")+
  theme_light()+
  theme(strip.text = element_text(size=20),plot.title = element_text(size=18, face="bold", vjust=2),axis.title.x = element_text(size=20,color="black", vjust=-0.5),
        axis.title.y = element_text(size=20,color="black" , vjust=0.5) )
g

## Function for arranging ggplots. use png(); arrange(p1, p2, ncol=1); dev.off() to save.
require(grid)
vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange_ggplot2 <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
  dots <- list(...)
  n <- length(dots)
  if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
  if(is.null(nrow)) { nrow = ceiling(n/ncol)}
  if(is.null(ncol)) { ncol = ceiling(n/nrow)}
  ## NOTE see n2mfrow in grDevices for possible alternative
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
  ii.p <- 1
  for(ii.row in seq(1, nrow)){
    ii.table.row <- ii.row	
    if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
    for(ii.col in seq(1, ncol)){
      ii.table <- ii.p
      if(ii.p > n) break
      print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
      ii.p <- ii.p + 1
    }
  }
}

save.image("fix.RData")
