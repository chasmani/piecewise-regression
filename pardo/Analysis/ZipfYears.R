
library(igraph)
library(MASS)
library(stats)
library(SOAR)
library(segmented)
library(psych)
library(ggplot2)
library(poweRlaw)
library(reshape2)
setwd("~/GoogleDrive/Project/Tables/")

#To load stored data
load(result, envir = parent.frame(), verbose = FALSE)
load("todo.RData")
varlist<-c("ger","spa","fre","ita","eng","rus","chi")
varlist_data<-lapply(varlist, function(c.name) {return(paste("data_", c.name, sep=""))})
#To remove all objects
rm(list=ls())
mapply(assign,varlist_data,lapply(varlist, function(c.name) {
  return(read.csv(paste("all_words_per_year_", c.name, ".txt", sep="")))}))

cit <- read.csv("citations.CSV", quote = "", 
                row.names = NULL, 
                stringsAsFactors = FALSE)
data_ger = read.csv("all_words_per_year_ger.txt", quote = "", row.names = NULL, stringsAsFactors = FALSE)
data_spa = read.csv("all_words_per_year_spa.txt", quote = "", row.names = NULL, stringsAsFactors = FALSE)
data_fre = read.csv("all_words_per_year_fre.txt", quote = "", row.names = NULL, stringsAsFactors = FALSE)
data_ita = read.csv("all_words_per_year_ita.txt", quote = "", row.names = NULL, stringsAsFactors = FALSE)
data_eng = read.csv("all_words_per_year_eng.txt", quote = "", row.names = NULL, stringsAsFactors = FALSE)
data_rus = read.csv("all_words_per_year_rus.txt", quote = "", row.names = NULL, stringsAsFactors = FALSE)



isPercentage<-function(var,total){
  if(var<2){
    var<-floor(total*var)
  }
  return(var)
}
calculateCoef<-function(freqs,init,end)
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
  total<-length(occur_per_year)
  #occur_per_year<-occur_per_year/total
  #The sequence of ranks is created
  x <- log(seq_along(occur_per_year))
  #The coefficients of the linear regressions are calculated
  if(length(occur_per_year)!=0)
  {  
    m <- lm(log(occur_per_year) ~ x)
    value<-coef(m)[[2]]
  }
  else
    value<-NA
  return(c(-value,total))
}

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
  occur_per_year<-occur_per_year/total
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
    #values1<-NA
    test<-NA
    test2<-NA
    xmin<-NA
  }
  return(c(-value,total,test,test2,xmin))

}

calculateCoefPowFix<-function(occur_per_year)
{
  
  #The words are ranked
  occur_per_year<-sort(occur_per_year, decreasing = TRUE)
  
  #Normalization
  total<-sum(occur_per_year)
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
calculateCoefFix<-function(data)
{
  totalData<-sum(as.numeric(data))
  maxDoc<-max(as.numeric(data))
  data<-data[data!=0]
  data<-data/totalData
  numWords<-length(data)
  data_hist<-try(hist(log(data),breaks = 30, plot=FALSE))
  if(class(data_hist) == "try-error") {return(c(NA,NA,NA,NA,NA))}
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

fix<-function(data)
{
  data<-sort(data,decreasing=TRUE)
  return(data[1:5000])
}
###Abreviatted Results
data<-list(data_ger,data_spa,data_fre,data_ita,data_eng,data_rus)#,data_chi)
dataFixKer<-lapply(data,function(name){return(apply(name[2:209],2,fix))})
names(dataFixKer)<-c("ger","spa","fre","ita","eng","rus")
resultsFix<-lapply(dataFixKer,function(name){return(apply(name[2:209],2,calculateCoefPowFix))})
names(resultsFix)<-c("ger","spa","fre","ita","eng","rus")
resultsNorm<-lapply(data,function(name){return(apply(name[2:209],2,calculateCoefPow,init=2,end=nrow(name)-1))})
names(resultsNorm)<-c("ger","spa","fre","ita","eng","rus")
total<-lapply(results,function(name){return(name[3,])})
names(total)<-c("ger","spa","fre","ita","eng","rus")
coefs<-lapply(results,function(name){return(name[1,])})
names(coefs)<-c("ger","spa","fre","ita","eng","rus")
coefsNorm<-lapply(results,function(name){return(name[2,])})
names(coefsNorm)<-c("ger","spa","fre","ita","eng","rus")
xmin<-lapply(results,function(name){return(name[5,])})
names(xmin)<-c("ger","spa","fre","ita","eng","rus")

totalFix<-lapply(resultsFix,function(name){return(name[2,])})
names(totalFix)<-c("ger","spa","fre","ita","eng","rus")
coefsFix<-lapply(resultsFix,function(name){return(name[1,])})
names(coefsFix)<-c("ger","spa","fre","ita","eng","rus")
testFix<-lapply(resultsFix,function(name){return(name[3,])})
names(testFix)<-c("ger","spa","fre","ita","eng","rus")
xminFix<-lapply(resultsFix,function(name){return(name[5,])})
names(xminFix)<-c("ger","spa","fre","ita","eng","rus")

segReg<-lapply(data,function(name){return(apply(name[2:209],2,segmentedRegression))})

#Plotting coefficients
year<-seq(1800,2007,1)
plot(year,coefsFix$ger,col="black",type="l",xlab="Year",ylab="Zipf's coefficient",main="Evolution of Zipf's coefficient for kernel lexicon")
points(year,coefsFix$spa,col="red",type="l")
points(year,coefsFix$fre,col="blue",type="l")
points(year,coefsFix$ita,col="green",type="l")
points(year,coefsFix$eng,col="purple",type="l")
points(year,coefsFix$rus,col="brown",type="l")
#points(year,coefsKer_chi,col="yellow",type="l")
legend("topleft", c("German","Spanish","French","Italian","English","Russian"),lty=c(1,1),col=c("black","red","blue","green","purple","brown","yellow"))

plot(year,coefsNorm$ger,col="black",type="l",xlab="Year",ylab="Zipf's coefficient",main="Evolution of Zipf's coefficient for kernel lexicon")
points(year,coefsNorm$spa,col="red",type="l")
points(year,coefsNorm$fre,col="blue",type="l")
points(year,coefsNorm$ita,col="green",type="l")
points(year,coefsNorm$eng,col="purple",type="l")
points(year,coefsNorm$rus,col="brown",type="l")
#points(year,coefsKer_chi,col="yellow",type="l")
legend("topleft", c("German","Spanish","French","Italian","English","Russian"),lty=c(1,1),col=c("black","red","blue","green","purple","brown","yellow"))

plot(year,xminFix$ger,col="black",type="l",xlab="Year",ylab="Zipf's coefficient",main="Evolution of Zipf's coefficient for kernel lexicon")
points(year,xminFix$spa,col="red",type="l")
points(year,xminFix$fre,col="blue",type="l")
points(year,xminFix$ita,col="green",type="l")
points(year,xminFix$eng,col="purple",type="l")
points(year,xminFix$rus,col="brown",type="l")
#points(year,coefsKer_chi,col="yellow",type="l")
legend("topleft", c("German","Spanish","French","Italian","English","Russian"),lty=c(1,1),col=c("black","red","blue","green","purple","brown","yellow"))

#
findKerLim<-function(data)
{
  xmin<-data[length(data)]
  data<-data[1:(length(data)-1)]
  data<-data[data!=0]
  return(length(data[data>xmin])/length(data))
}
findKerLimAbs<-function(data)
{
  xmin<-data[length(data)]
  data<-data[1:(length(data)-1)]
  data<-data[data!=0]
  return(length(data[data>xmin]))
}

dataSpa<-rbind(data_spa[,2:209],xmin$spa)
kerLimSpa<-apply(dataSpa,2,findKerLim)
plot(year,kerLimSpa,xlab="Year",ylab="xmin",main="Evolution of relative xmin for spanish")
plot(year,xmin$spa,xlab="Year",ylab="xmin",main="Evolution of absolute xmin for spanish")

dataSpaFix<-rbind(dataFixKer$spa,xminFix$spa)
kerLimSpaFix<-apply(dataSpaFix,2,findKerLim)
plot(year,kerLimSpaFix,xlab="Year",ylab="xmin",main="Evolution of relative xmin for spanish")
plot(year,xmin$spa,xlab="Year",ylab="xmin",main="Evolution of absolute xmin for spanish")

dataGer<-rbind(data_ger[,2:209],xmin$ger)
kerLimGer<-apply(dataGer,2,findKerLim)
plot(year,kerLimGer,xlab="Year",ylab="xmin",main="Evolution of xmin for german")
plot(year,xmin$ger,xlab="Year",ylab="xmin",main="Evolution of absolute xmin for german")

dataSpaFix<-rbind(dataFixKer$spa,xminFix$spa)
kerLimSpaFix<-apply(dataSpaFix,2,findKerLim)
plot(year,kerLimSpaFix,xlab="Year",ylab="xmin",main="Evolution of relative xmin for spanish")
plot(year,xmin$spa,xlab="Year",ylab="xmin",main="Evolution of absolute xmin for spanish")

dataFre<-rbind(data_fre[,2:209],xmin$fre)
kerLimFre<-apply(dataFre,2,findKerLim)
plot(year,kerLimFre,xlab="Year",ylab="xmin",main="Evolution of xmin for french")
plot(year,xmin$fre,xlab="Year",ylab="xmin",main="Evolution of absolute xmin for french")

dataFreFix<-rbind(dataFixKer$fre,xminFix$fre)
kerLimFreFix<-apply(dataFreFix,2,findKerLim)
plot(year,kerLimFreFix,xlab="Year",ylab="xmin",main="Evolution of relative xmin for spanish")
plot(year,xmin$spa,xlab="Year",ylab="xmin",main="Evolution of absolute xmin for spanish")

dataIta<-rbind(data_ita[,2:209],xmin$ita)
kerLimIta<-apply(dataIta,2,findKerLim)
plot(year,kerLimIta,xlab="Year",ylab="xmin",main="Evolution of xmin for italian")
plot(year,xmin$ita,xlab="Year",ylab="xmin",main="Evolution of absolute xmin for italian")

dataEng<-rbind(data_eng[,2:209],xmin$eng)
kerLimEng<-apply(dataEng,2,findKerLim)
plot(year,kerLimEng,xlab="Year",ylab="xmin",main="Evolution of xmin for english")
plot(year,xmin$ita,xlab="Year",ylab="xmin",main="Evolution of absolute xmin for italian")

dataEngFix<-rbind(dataFixKer$eng,xminFix$eng)
kerLimEngFix<-apply(dataEngFix,2,findKerLim)
plot(year,kerLimEngFix,xlab="Year",ylab="xmin",main="Evolution of relative xmin for spanish")
plot(year,xmin$spa,xlab="Year",ylab="xmin",main="Evolution of absolute xmin for spanish")

dataRus<-rbind(data_rus[,2:209],xmin$rus)
kerLimRus<-apply(dataRus,2,findKerLim)
plot(year,kerLimRus,xlab="Year",ylab="xmin",main="Evolution of xmin for russian")
plot(year,xmin$ita,xlab="Year",ylab="xmin",main="Evolution of absolute xmin for italian")

totalts<-lapply(results,function(name){return(ts(name[2,][!is.na(name[1,])&!is.na(name[2,])]))})
coefsts<-lapply(results,function(name){return(ts(name[1,][!is.na(name[1,])&!is.na(name[2,])]))})
xmin<-lapply(results,function(name){return(ts(name[5,][!is.na(name[1,])&!is.na(name[2,])]))})


year<-seq(1800,2007,1)
varlist_coefs<-lapply(varlist, function(c.name) {return(paste("coefs_", c.name, sep=""))})
varlist_coefs<-lapply(varlist,function(c.name){})
lapply(varlist1, function(c.name, x) { 
    png(paste("c:/TimePlot-", c.name, ".png", sep=""))
    plot(x[,c.name], main=c.name)
    dev.off()
  }, x=x)




##All data
results_ger<-apply(data_ger[2:211],2,calculateCoef,init=1,end=nrow(data_ger)-1)
coefs_ger<-results_ger[1,]
total_ger<-results_ger[2,]
results_spa<-apply(data_spa[2:211],2,calculateCoef,init=1,end=nrow(data_spa)-1)
coefs_spa<-results_spa[1,]
total_spa<-results_spa[2,]
results_fre<-apply(data_fre[2:211],2,calculateCoef,init=1,end=nrow(data_fre)-1)
coefs_fre<-results_fre[1,]
total_fre<-results_fre[2,]
results_ita<-apply(data_ita[2:211],2,calculateCoef,init=1,end=nrow(data_ita)-1)
coefs_ita<-results_ita[1,]
total_ita<-results_ita[2,]
results_eng<-apply(data_eng[3:201],2,calculateCoef,init=1,end=nrow(data_eng)-1)
coefs_eng<-results_eng[1,]
total_eng<-results_eng[2,]
coefs_engtemp<-rep(0,210)
coefs_engtemp[2:200]<-coefs_eng
coefs_engtemp[1]<-NA
coefs_engtemp[201:210]<-NA
coefs_eng<-coefs_engtemp

total_engtemp<-rep(0,210)
total_engtemp[2:200]<-total_eng
total_engtemp[1]<-NA
total_engtemp[201:210]<-NA
total_eng<-total_engtemp

results_rus<-apply(data_rus[2:211],2,calculateCoef,init=1,end=nrow(data_rus)-1)
coefs_rus<-results_rus[1,]
total_rus<-results_rus[2,]

resultsPow_eng<-apply(data_eng[2:209],2,calculateCoefPow,init=2,end=nrow(data_eng)-1)
coefsPow_eng<-resultsPow_eng[1,]
totalPow_eng<-resultsPow_eng[2,]
testPow_eng<-resultsPow_eng[3,]
test2Pow_eng<-resultsPow_eng[4,]
xminPow_eng<-resultsPow_eng[5,]

resultsPow_fre<-apply(data_fre[2:209],2,calculateCoefPow,init=2,end=nrow(data_fre)-1)
coefsPow_fre<-resultsPow_fre[1,]
totalPow_fre<-resultsPow_fre[2,]
testPow_fre<-resultsPow_fre[3,]
test2Pow_fre<-resultsPow_fre[4,]
xminPow_fre<-resultsPow_fre[5,]

resultsPow_spa<-apply(data_spa[2:209],2,calculateCoefPow,init=2,end=nrow(data_spa)-1)
coefsPow_spa<-resultsPow_spa[1,]
totalPow_spa<-resultsPow_spa[2,]
testPow_spa<-resultsPow_spa[3,]
test2Pow_spa<-resultsPow_spa[4,]
xminPow_spa<-resultsPow_spa[5,]

resultsPow_ger<-apply(data_ger[2:209],2,calculateCoefPow,init=2,end=nrow(data_ger)-1)
coefsPow_ger<-resultsPow_ger[1,]
totalPow_ger<-resultsPow_ger[2,]
testPow_ger<-resultsPow_ger[3,]
test2Pow_ger<-resultsPow_ger[4,]
xminPow_ger<-resultsPow_ger[5,]

results_chi<-apply(data_chi[2:211],2,calculateCoef,init=1,end=nrow(data_chi)-1)
coefs_chi<-results_chi[1,]
total_chi<-results_chi[2,]
###Time Series
total_gerts<-ts(total$ger[!is.na(coefs$ger)&!is.na(total$ger)])
total_spats<-ts(total$spa[!is.na(coefs$spa)&!is.na(total$spa)])
total_frets<-ts(total$fre[!is.na(coefs$fre)&!is.na(total$fre)])
total_itats<-ts(total$ita[!is.na(coefs$ita)&!is.na(total$ita)])
total_engts<-ts(total$eng[!is.na(coefs$eng)&!is.na(total$eng)])
total_rusts<-ts(total$rus[!is.na(coefs$rus)&!is.na(total$rus)])

xmin_gerts<-ts(xmin$ger[!is.na(coefs$ger)&!is.na(xmin$ger)])
xmin_spats<-ts(xmin$spa[!is.na(coefs$spa)&!is.na(xmin$spa)])
xmin_frets<-ts(xmin$fre[!is.na(coefs$fre)&!is.na(xmin$fre)])
xmin_itats<-ts(xmin$ita[!is.na(coefs$ita)&!is.na(xmin$ita)])
xmin_engts<-ts(xmin$eng[!is.na(coefs$eng)&!is.na(xmin$eng)])
xmin_rusts<-ts(xmin$rus[!is.na(coefs$rus)&!is.na(xmin$rus)])

##All data

coefs_gerts<-ts(coefs$ger[!is.na(coefs$ger)&!is.na(xmin$ger)])
coefs_spats<-ts(coefs$spa[!is.na(coefs$spa)&!is.na(xmin$spa)])
coefs_frets<-ts(coefs$fre[!is.na(coefs$fre)&!is.na(xmin$fre)])
coefs_itats<-ts(coefs$ita[!is.na(coefs$ita)&!is.na(xmin$ita)])
coefs_engts<-ts(coefs$eng[!is.na(coefs$eng)&!is.na(xmin$eng)])
coefs_rusts<-ts(coefs$rus[!is.na(coefs$rus)&!is.na(xmin$rus)])

ccfvalues<-mapply(function(name1,name2){return(ccf(name1,name2,0,plot=FALSE))},coefs,total)
ccfvalues_ger = ccf(xmin_gerts,coefs_gerts,0,plot=FALSE)
ccfvalues_ger
ccfvalues_spa = ccf(xmin_spats,coefs_spats,0,plot=FALSE)
ccfvalues_spa
ccfvalues_fre = ccf(xmin_frets,coefs_frets,0,plot=FALSE)
ccfvalues_fre
ccfvalues_ita = ccf(xmin_itats,coefs_itats,0,plot=FALSE)
ccfvalues_ita
ccfvalues_eng = ccf(xmin_engts,coefs_engts,0,plot=FALSE)
ccfvalues_eng
ccfvalues_rus = ccf(xmin_rusts,coefs_rusts,0,plot=FALSE)
ccfvalues_rus

##Plot
year<-seq(1850,2007,1)
plot(year,xminPow_ger[51:length(xminPow_ger)],col="black",type="l",xlab="Year",ylab="Zipf's coefficient",main="Evolution of Zipf's coefficient")
points(year,xminPow_spa[51:length(xminPow_spa)],col="red",type="l")
points(year,xminPow_fre[51:length(xminPow_fre)],col="blue",type="l")
points(year,coefs_ita,col="green",type="l")
points(year,coefs_eng,col="purple",type="l")
points(year,coefsPow_rus,col="brown",type="l")
points(year,coefs_chi,col="yellow",type="l")
legend("topleft", c("German","Spanish","French","Italian","English","Russian","Chinese"),lty=c(1,1),col=c("black","red","blue","green","purple","brown","yellow"))
m_ger<-lm(coefs_ger~year)
m_spa<-lm(coefs_spa~year)
m_fre<-lm(coefs_fre~year)
m_ita<-lm(coefs_ita~year)
abline(m_ger,col="black")
abline(m_spa,col="red")
abline(m_fre,col="blue")
abline(m_ita,col="green")

#Total vs Time
plot(year,total_rus,col="purple",type="l",xlab="Year",ylab="Zipf's coefficient",main="Evolution of Number of Different Words")
points(year,total_spa,col="red",type="l")
points(year,total_fre,col="blue",type="l")
points(year,total_ita,col="green",type="l")
points(year,total_ger,col="black",type="l")
points(year,total_rus,col="brown",type="l")
points(year,total_chi,col="yellow",type="l")
legend("topleft", c("German","Spanish","French","Italian","English","Russian","Chinese"),lty=c(1,1),col=c("black","red","blue","green","purple","brown","yellow"),cex = 0.75)

#Total number of words vs Zipf

plot(total_ger,coefs_ger,col="black",type="l",xlab="Year",ylab="Zipf's coefficient",ylim=c(1.1,1.8),,main="Evolution of Zipf's coefficient")
points(total_spa,coefs_spa,col="red",type="l")
points(total_fre,coefs_fre,col="blue",type="l")
points(total_ita,coefs_ita,col="green",type="l")
legend("topleft", c("German","Spanish","French","Italian"),lty=c(1,1),col=c("black","red","blue","green"))
m_ger<-lm(coefs_ger~total_ger)
m_spa<-lm(coefs_spa~total_spa)
m_fre<-lm(coefs_fre~total_fre)
m_ita<-lm(coefs_ita~total_ita)
abline(m_ger,col="black")
abline(m_spa,col="red")
abline(m_fre,col="blue")
abline(m_ita,col="green")

list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}


##Popular words


results_ger<-apply(data_ger[2:211],2,calculateCoef,init=1,end=0.05)
coefsKer_ger<-results_ger[1,]
totalKer_ger<-results_ger[2,]
results_spa<-apply(data_spa[2:211],2,calculateCoef,init=1,end=0.05)
coefsKer_spa<-results_spa[1,]
totalKer_spa<-results_spa[2,]
results_fre<-apply(data_fre[2:211],2,calculateCoef,init=1,end=0.05)
coefsKer_fre<-results_fre[1,]
totalKer_fre<-results_fre[2,]
results_ita<-apply(data_ita[2:211],2,calculateCoef,init=1,end=0.05)
coefsKer_ita<-results_ita[1,]
totalKer_ita<-results_ita[2,]
results_eng<-apply(data_eng[3:201],2,calculateCoef,init=1,end=0.05)
coefsKer_eng<-results_eng[1,]
totalKer_eng<-results_eng[2,]
coefsKer_engtemp<-rep(0,210)
coefsKer_engtemp[2:200]<-coefsKer_eng
coefsKer_engtemp[1]<-NA
coefsKer_engtemp[201:210]<-NA
coefsKer_eng<-coefsKer_engtemp
totalKer_engtemp<-rep(0,210)
totalKer_engtemp[2:200]<-totalKer_eng
totalKer_engtemp[1]<-NA
totalKer_engtemp[201:210]<-NA
totalKer_eng<-totalKer_engtemp

results_rus<-apply(data_rus[2:211],2,calculateCoef,init=1,end=0.05)
coefsKer_rus<-results_rus[1,]
totalKer_rus<-results_rus[2,]
results_chi<-apply(data_chi[2:211],2,calculateCoef,init=1,end=0.05)
coefsKer_chi<-results_chi[1,]
totalKer_chi<-results_chi[2,]

plot(year,coefsKer_ger,col="black",type="l",xlab="Year",ylab="Zipf's coefficient",ylim=c(0.7,1.4),main="Evolution of Zipf's coefficient for kernel lexicon")
points(year,coefsKer_spa,col="red",type="l")
points(year,coefsKer_fre,col="blue",type="l")
points(year,coefsKer_ita,col="green",type="l")
points(year,coefsKer_eng,col="purple",type="l")
points(year,coefsKer_rus,col="brown",type="l")
points(year,coefsKer_chi,col="yellow",type="l")
legend("topleft", c("German","Spanish","French","Italian","English","Russian","Chinese"),lty=c(1,1),col=c("black","red","blue","green","purple","brown","yellow"))
m_ger<-lm(coefsKer_ger~year)
m_spa<-lm(coefsKer_spa~year)
m_fre<-lm(coefsKer_fre~year)
m_ita<-lm(coefsKer_ita~year)
abline(m_ger,col="black")
abline(m_spa,col="red")
abline(m_fre,col="blue")
abline(m_ita,col="green")

###Time Series
totalKer_gerts<-ts(totalKer_ger[!is.na(coefsKer_ger)&!is.na(totalKer_ger)])
totalKer_spats<-ts(totalKer_spa[!is.na(coefsKer_spa)&!is.na(totalKer_spa)])
totalKer_frets<-ts(totalKer_fre[!is.na(coefsKer_fre)&!is.na(totalKer_fre)])
totalKer_itats<-ts(totalKer_ita[!is.na(coefsKer_ita)&!is.na(totalKer_ita)])
totalKer_engts<-ts(totalKer_eng[!is.na(coefsKer_eng)&!is.na(totalKer_eng)])
totalKer_rusts<-ts(totalKer_rus[!is.na(coefsKer_rus)&!is.na(totalKer_rus)])
totalKer_chits<-ts(totalKer_chi[!is.na(coefsKer_chi)&!is.na(totalKer_chi)])

##All data
coefsKer_gerts<-ts(coefsKer_ger[!is.na(coefsKer_ger)&!is.na(totalKer_ger)])
coefsKer_spats<-ts(coefsKer_spa[!is.na(coefsKer_spa)&!is.na(totalKer_spa)])
coefsKer_frets<-ts(coefsKer_fre[!is.na(coefsKer_fre)&!is.na(totalKer_fre)])
coefsKer_itats<-ts(coefsKer_ita[!is.na(coefsKer_ita)&!is.na(totalKer_ita)])
coefsKer_engts<-ts(coefsKer_eng[!is.na(coefsKer_eng)&!is.na(totalKer_eng)])
coefsKer_rusts<-ts(coefsKer_rus[!is.na(coefsKer_rus)&!is.na(totalKer_rus)])
coefsKer_chits<-ts(coefsKer_chi[!is.na(coefsKer_chi)&!is.na(totalKer_chi)])

ccfvaluesKer_ger = ccf(totalKer_gerts,coefsKer_gerts,0,plot=FALSE)
ccfvaluesKer_ger
ccfvaluesKer_spa = ccf(totalKer_spats,coefsKer_spats,0,plot=FALSE)
ccfvaluesKer_spa
ccfvaluesKer_fre = ccf(totalKer_frets,coefsKer_frets,0,plot=FALSE)
ccfvaluesKer_fre
ccfvaluesKer_ita = ccf(totalKer_itats,coefsKer_itats,0,plot=FALSE)
ccfvaluesKer_ita
ccfvaluesKer_eng = ccf(totalKer_engts,coefsKer_engts,0,plot=FALSE)
ccfvaluesKer_eng
ccfvaluesKer_rus = ccf(totalKer_rusts,coefsKer_rusts,0,plot=FALSE)
ccfvaluesKer_rus
ccfvaluesKer_chi = ccf(totalKer_chits,coefsKer_chits,0,plot=FALSE)
ccfvaluesKer_chi




#Power law
coefsKerPow<-apply(mydata[2:201],2,calculateCoefPow,init=1,end=500)
plot(x,coefsKerPow)
m<-lm(coefsKerPow~x)
abline(m)



##Unlimited lexicon
#LR
##All data
results_ger<-apply(data_ger[2:211],2,calculateCoef,init=0.05,end=nrow(data_ger)-1)
coefsUnl_ger<-results_ger[1,]
totalUnl_ger<-results_ger[2,]
results_spa<-apply(data_spa[2:211],2,calculateCoef,init=0.05,end=nrow(data_spa)-1)
coefsUnl_spa<-results_spa[1,]
totalUnl_spa<-results_spa[2,]
results_fre<-apply(data_fre[2:211],2,calculateCoef,init=0.05,end=nrow(data_fre)-1)
coefsUnl_fre<-results_fre[1,]
totalUnl_fre<-results_fre[2,]
results_ita<-apply(data_ita[2:211],2,calculateCoef,init=0.05,end=nrow(data_ita)-1)
coefsUnl_ita<-results_ita[1,]
totalUnl_ita<-results_ita[2,]
results_eng<-apply(data_eng[3:201],2,calculateCoef,init=0.05,end=nrow(data_eng)-1)
coefsUnl_eng<-results_eng[1,]
totalUnl_eng<-results_eng[2,]
coefsUnl_engtemp<-rep(0,210)
coefsUnl_engtemp[2:200]<-coefsUnl_eng
coefsUnl_engtemp[1]<-NA
coefsUnl_engtemp[201:210]<-NA
coefsUnl_eng<-coefsUnl_engtemp

totalUnl_engtemp<-rep(0,210)
totalUnl_engtemp[2:200]<-total_eng
totalUnl_engtemp[1]<-NA
totalUnl_engtemp[201:210]<-NA
totalUnl_eng<-total_engtemp

results_rus<-apply(data_rus[2:211],2,calculateCoef,init=0.05,end=nrow(data_rus)-1)
coefsUnl_rus<-results_rus[1,]
totalUnl_rus<-results_rus[2,]
results_chi<-apply(data_chi[2:211],2,calculateCoef,init=0.05,end=nrow(data_chi)-1)
coefsUnl_chi<-results_chi[1,]
totalUnl_chi<-results_chi[2,]
###Time Series
totalUnl_gerts<-ts(totalUnl_ger[!is.na(coefsUnl_ger)&!is.na(totalUnl_ger)])
totalUnl_spats<-ts(totalUnl_spa[!is.na(coefsUnl_spa)&!is.na(totalUnl_spa)])
totalUnl_frets<-ts(totalUnl_fre[!is.na(coefsUnl_fre)&!is.na(totalUnl_fre)])
totalUnl_itats<-ts(totalUnl_ita[!is.na(coefsUnl_ita)&!is.na(totalUnl_ita)])
totalUnl_engts<-ts(totalUnl_eng[!is.na(coefsUnl_eng)&!is.na(totalUnl_eng)])
totalUnl_rusts<-ts(totalUnl_rus[!is.na(coefsUnl_rus)&!is.na(totalUnl_rus)])
totalUnl_chits<-ts(totalUnl_chi[!is.na(coefsUnl_chi)&!is.na(totalUnl_chi)])

##All data
coefsUnl_gerts<-ts(coefsUnl_ger[!is.na(coefsUnl_ger)&!is.na(totalUnl_ger)])
coefsUnl_spats<-ts(coefsUnl_spa[!is.na(coefsUnl_spa)&!is.na(totalUnl_spa)])
coefsUnl_frets<-ts(coefsUnl_fre[!is.na(coefsUnl_fre)&!is.na(totalUnl_fre)])
coefsUnl_itats<-ts(coefsUnl_ita[!is.na(coefsUnl_ita)&!is.na(totalUnl_ita)])
coefsUnl_engts<-ts(coefsUnl_eng[!is.na(coefsUnl_eng)&!is.na(totalUnl_eng)])
coefsUnl_rusts<-ts(coefsUnl_rus[!is.na(coefsUnl_rus)&!is.na(totalUnl_rus)])
coefsUnl_chits<-ts(coefsUnl_chi[!is.na(coefsUnl_chi)&!is.na(totalUnl_chi)])

ccfvaluesUnl_ger = ccf(totalUnl_gerts,coefsUnl_gerts,0,plot=FALSE)
ccfvaluesUnl_ger
ccfvaluesUnl_spa = ccf(totalUnl_spats,coefsUnl_spats,0,plot=FALSE)
ccfvaluesUnl_spa
ccfvaluesUnl_fre = ccf(totalUnl_frets,coefsUnl_frets,0,plot=FALSE)
ccfvaluesUnl_fre
ccfvaluesUnl_ita = ccf(totalUnl_itats,coefsUnl_itats,0,plot=FALSE)
ccfvaluesUnl_ita
ccfvaluesUnl_eng = ccf(totalUnl_engts,coefsUnl_engts,0,plot=FALSE)
ccfvaluesUnl_eng
ccfvaluesUnl_rus = ccf(totalUnl_rusts,coefsUnl_rusts,0,plot=FALSE)
ccfvaluesUnl_rus
ccfvaluesUnl_chi = ccf(totalUnl_chits,coefsUnl_chits,0,plot=FALSE)
ccfvaluesUnl_chi

##Plot
plot(year,coefsUnl_ger,col="black",type="l",xlab="Year",ylab="Zipf's coefficient",ylim=c(0.5,1.90),main="Evolution of Zipf's coefficient for the unlimited lexicon")
points(year,coefsUnl_spa,col="red",type="l")
points(year,coefsUnl_fre,col="blue",type="l")
points(year,coefsUnl_ita,col="green",type="l")
points(year,coefsUnl_eng,col="purple",type="l")
points(year,coefsUnl_rus,col="brown",type="l")
points(year,coefsUnl_chi,col="yellow",type="l")
legend("topleft", c("German","Spanish","French","Italian","English","Russian","Chinese"),lty=c(1,1),col=c("black","red","blue","green","purple","brown","yellow"))
m_ger<-lm(coefsUnl_ger~year)
m_spa<-lm(coefsUnl_spa~year)
m_fre<-lm(coefsUnl_fre~year)
m_ita<-lm(coefsUnl_ita~year)
abline(m_ger,col="black")
abline(m_spa,col="red")
abline(m_fre,col="blue")
abline(m_ita,col="green")
#PL
coefsUnlPow<-apply(mydata[2:201],2,calculateCoefPow,init=501,end=nrow(mydata)-1)
plot(x,coefsUnlPow)
m<-lm(coefsUnlPow~x)
abline(m)



####Plotting one year

year1<-data_fre[200][!is.na(data_fre[200])]
year1<-year1[year1!=0]
year1<-sort(year1, decreasing = TRUE)
x <- log(seq_along(year1))
powfit<-power.law.fit(year1,implementation='plfit')
yearLim<-year1[year1>powfit$xmin]
xLim<-log(seq_along(yearLim))
plot(x,log(year1))
plot(xLim,log(yearLim))
abline(14,-powfit$alpha)
length(yearLim)/length(year1)

year2<-data_spa[100][!is.na(data_spa[100])]
year2<-year2[year2!=0]
year2<-sort(year2, decreasing = TRUE)
x <- log(seq_along(year2))
powfit<-power.law.fit(year2,implementation='plfit')
plot(x,log(year2))
abline(14,-powfit$alpha)

data<-data_fre[,52:209]
data<-rbind(data,xminPow_fre[51:208])

findKerLim<-function(data)
{
  xmin<-data[length(data)]
  data<-data[1:(length(data)-1)]
  data<-data[data!=0]
  return(length(data[data>xmin])/length(data))
}

kerLim<-apply(data,2,findKerLim)
year<-seq(1850,2007,1)
plot(year,kerLim)


data1<-data_spa[200]
totalData<-sum(data1)
maxDoc<-max(data1)
data1<-data1[data1!=0]
data1<-data1/totalData
numWords<-length(data1)
data_hist<-hist(log(data1),probability = FALSE, breaks = 20, col = "darkslategray4", border = "seashell3")
x<-sort(data_hist$counts[data_hist$counts!=0],decreasing=FALSE)
y<-sort(exp(data_hist$mids[data_hist$counts!=0]),decreasing=TRUE)
x<-x[1:(length(x)-2)]
y<-y[1:(length(y)-2)]
#Plot empirical vs fitted distribution density functions
x<-log(x)
y<-log(y)
lin.mod<-lm(y~x)
segmented.mod <- segmented(lin.mod, seg.Z = ~x, psi=7)
plot(x,y)
plot(segmented.mod, add=T)
summary(segmented.mod)
a<-power.law.fit(x,implementation="plfit")
x1<-sort(data_hist$counts[data_hist$counts<=a$xmin],decreasing=FALSE)
y1<-sort(exp(data_hist$mids[data_hist$counts<=a$xmin]),decreasing=TRUE)
plot(x1,y1,log="xy")
#recovering the break point
exp(segmented.mod$psi[2])

a<-data_hist$counts[data_hist$counts>exp(segmented.mod$psi[2])]
b<-data_hist$mids[data_hist$counts>exp(segmented.mod$psi[2])]
b<-exp(b)
b<-b*totalData
c<-a*b
sum(c)/totalData

segmentedRegression<-function(data)
{
  totalData<-sum(data)
  maxDoc<-max(data)
  data<-data[data!=0]
  data<-data/totalData
  numWords<-length(data)
  data_hist<-try(hist(log(data),breaks = 20, plot=FALSE))
  if(class(data_hist) == "try-error") {return(c(NA,NA,NA,NA,NA))}
  x<-sort(data_hist$counts[data_hist$counts!=0],decreasing=FALSE)
  y<-sort(exp(data_hist$mids[data_hist$counts!=0]),decreasing=TRUE)
  x<-x[1:(length(x)-2)]
  y<-y[1:(length(y)-2)]
  #Where the power law stops holding
  aux<-power.law.fit(y,implementation="plfit")$xmin
  #Plot empirical vs fitted distribution density functions
  x<-log(x)
  y<-log(y)
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
segReg<-lapply(data,function(name){return(apply(name[2:209],2,segmentedRegression))})
data1<-data_ger[20]
reg<-segmentedRegression(data1) #apply(data1,2,segmentedRegression,init=2,end=nrow(name)-1))
reg
d4<-apply(data_fre[2:209],2,segmentedRegression)
year<-seq(1800,2007,1)
plot(year,d4[1,],col="black",type="l",xlab="Year",ylab="log(xmin)",main="Evolution of $log(xmin)$")
plot(year,d4[2,],col="black",type="l",xlab="Year",ylab="log(xmin)",main="Evolution of $log(xmin)$")
plot(year,d4[3,],col="black",type="l",xlab="Year",ylab="log(xmin)",main="Evolution of $log(xmin)$")
plot(year,d4[4,],col="black",type="l",xlab="Year",ylab="log(xmin)",main="Evolution of $log(xmin)$")
plot(year,d4[5,],col="black",type="l",xlab="Year",ylab="log(xmin)",main="Evolution of $log(xmin)$")
plot(year,d4[6,],col="black",type="l",xlab="Year",ylab="log(xmin)",main="Evolution of $log(xmin)$")

#spa->4,10-30?
#ger->18
#ita->Somewhere before 20 but with histogram,above 30
names(segReg)<-c("ger","spa","fre","ita")#,"eng","rus","chi")
year<-seq(1848,2007,1)
plot(year,segReg$ger[1,],col="black",type="l",xlab="Year",ylab="Relative size of kernel lexicon",main="Evolution of the relative size of kernel lexicon")
plot(year,segReg$spa[1,],col="red",type="l")
plot(year,segReg$fre[1,],col="blue",type="l")
plot(year,segReg$ita[1,],col="green",type="l")
KL_ger<-lm(segReg$ger[1,]~year)
KL_spa<-lm(segReg$spa[1,]~year)
KL_fre<-lm(segReg$fre[1,]~year)
KL_ita<-lm(segReg$ita[1,]~year)
abline(KL_ger,col="black")
abline(KL_spa,col="red")
abline(KL_fre,col="blue")
abline(KL_ita,col="green")
plot(year,segReg$ger[2,],col="black",type="l",xlab="Year",ylab="alpha 1",main="Evolution of alpha 1")
plot(year,segReg$spa[2,],col="red",type="l")
plot(year,segReg$fre[2,],col="blue",type="l")
plot(year,segReg$ita[2,],col="green",type="l")
a1_ger<-lm(segReg$ger[2,]~year)
a1_spa<-lm(segReg$spa[2,]~year)
a1_fre<-lm(segReg$fre[2,]~year)
a1_ita<-lm(segReg$ita[2,]~year)
abline(a1_ger,col="black")
abline(a1_spa,col="red")
abline(a1_fre,col="blue")
abline(a1_ita,col="green")
plot(year,segReg$ger[3,],col="black",type="l",xlab="Year",ylab="alpha 2",main="Evolution of alpha 2")
points(year,segReg$spa[3,],col="red",type="l")
points(year,segReg$fre[3,],col="blue",type="l")
points(year,segReg$ita[3,],col="green",type="l")
a2_ger<-lm(segReg$ger[3,]~year)
a2_spa<-lm(segReg$spa[3,]~year)
a2_fre<-lm(segReg$fre[3,]~year)
a2_ita<-lm(segReg$ita[3,]~year)
abline(a2_ger,col="black")
abline(a2_spa,col="red")
abline(a2_fre,col="blue")
abline(a2_ita,col="green")
plot(year,segReg$ger[4,],col="black",type="l",xlab="Year",ylab="log(xmin)",main="Evolution of $log(xmin)$")
plot(year,segReg$spa[4,],col="red",type="l")
plot(year,segReg$fre[4,],col="blue",type="l")
plot(year,segReg$ita[4,],col="green",type="l")
xmin_ger<-lm(segReg$ger[4,]~year)
xmin_spa<-lm(segReg$spa[4,]~year)
xmin_fre<-lm(segReg$fre[4,]~year)
xmin_ita<-lm(segReg$ita[4,]~year)
abline(xmin_ger,col="black")
abline(xmin_spa,col="red")
abline(xmin_fre,col="blue")
abline(xmin_ita,col="green")
plot(year,segReg$ger[5,],col="black",type="l",xlab="Year",ylab="Kink",main="Evolution of the kink")
plot(year,segReg$spa[5,],col="red",type="l")
plot(year,segReg$fre[5,],col="blue",type="l")
plot(year,segReg$ita[5,],col="green",type="l")
kink_ger<-lm(segReg$ger[5,]~year)
kink_spa<-lm(segReg$spa[5,]~year)
kink_fre<-lm(segReg$fre[5,]~year)
kink_ita<-lm(segReg$ita[5,]~year)
abline(kink_ger,col="black")
abline(kink_spa,col="red")
abline(kink_fre,col="blue")
abline(kink_ita,col="green")
#points(year,coefsKer_eng,col="purple",type="l")
#points(year,coefsKer_rus,col="brown",type="l")
#points(year,coefsKer_chi,col="yellow",type="l")
legend("topleft", c("German","Spanish","French","Italian"),lty=c(1,1),col=c("black","red","blue","green"))#,"purple","brown","yellow"))

#Function that finds the index of the place in the cumulative sum of words to which the 
#random number corresponds
findPosition<-function(num,cumVector)
{
  return(which(cumVector>=num)[1])
}
#Function that performs a bootstrap over the data.
bootstrapData<-function(data)
{
  n<-1000000
  data<-data_eng[[i]]
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
#Function that calculates the quantities for one bootstrap of the data
statBootstrap<-function(data)
{
  newData<-replicate(100,bootstrapData(data))
  totalStats<-apply(newData,2,segmentedRegression)
  #Calculate the average of the bootstrap statistics
  #averages<-apply(totalStats,1,mean)
  #print("aaaaaa")
  return(totalStats)
}

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
g<-statBootstrap(data_ger[50])

finalResultsFre<-apply(data_fre[2:209],2,statBootstrap)
finalResultsSpa<-apply(data_spa[2:209],2,statBootstrap)
finalResultsGer<-apply(data_ger[2:209],2,statBootstrap)
finalResultsIta<-apply(data_ita[2:209],2,statBootstrap)
finalResultsEng<-apply(data_eng[2:209],2,statBootstrap)
finalResultsRus<-apply(data_rus[2:209],2,statBootstrap)

finalResultsFreFix<-apply(data_fre[2:209],2,statBootstrapFix)
finalResultsSpaFix<-apply(data_spa[2:209],2,statBootstrapFix)

separatingResults<-function(data,n)
{
  return(data[seq(1,600,1)%%6==n])
}
relativeKerFre<-as.data.frame(apply(finalResultsFre,2,separatingResults,n=1))
alpha1Fre<-as.data.frame(apply(finalResultsFre,2,separatingResults,n=2))
alpha2Fre<-as.data.frame(apply(finalResultsFre,2,separatingResults,n=3))
xminFre<-as.data.frame(apply(finalResultsFre,2,separatingResults,n=4))
kinkFre<-as.data.frame(apply(finalResultsFre,2,separatingResults,n=5))
psi.initFre<-as.data.frame(apply(finalResultsFre,2,separatingResults,n=0))

boxplot(relativeKerFre,outline=FALSE,xlab="year",ylab="% of vocabulary in the kernel lexicon",main=expression("Evolution of the relative size of the kernel lexicon for French"),names=as.character(seq(1800,2007,1)))
boxplot(alpha1Fre,outline=FALSE,xlab="year",ylim=c(0.95,1.2),ylab=expression(alpha[1]),main=expression(paste("Evolution of ", alpha[1], " for French")))
boxplot(alpha2Fre,outline=FALSE,ylim=c(1,3),xlab="year",ylab=expression(alpha[2]),main=expression(paste("Evolution of ",alpha[2]," for French")))
boxplot(xminFre,outline=FALSE,xlab="year",ylab="xmin",main=expression("Evolution of xmin for French"))
boxplot(kinkFre,outline=FALSE,xlab="year",ylab="Kink",main=expression("Evolution of the kink for French"))
boxplot(psi.initFre,outline=FALSE,xlab="year",ylab="Kink",main=expression("Evolution of the kink for French"))

boxplot(finalResultsFreFix,outline=FALSE,xlab="year",ylab=expression(alpha[1]),main=expression(paste("Evolution of ", alpha[1], " for French")))

relativeKerSpa<-as.data.frame(apply(finalResultsSpa,2,separatingResults,n=1))
alpha1Spa<-as.data.frame(apply(finalResultsSpa,2,separatingResults,n=2))
alpha2Spa<-as.data.frame(apply(finalResultsSpa,2,separatingResults,n=3))
xminSpa<-as.data.frame(apply(finalResultsSpa,2,separatingResults,n=4))
kinkSpa<-as.data.frame(apply(finalResultsSpa,2,separatingResults,n=0))
boxplot(relativeKerSpa,outline=FALSE,xlab="year",ylab="% of vocabulary in the kernel lexicon",main=expression("Evolution of the relative size of the kernel lexicon for Spanish"),names=as.character(seq(1800,2007,1)))
boxplot(alpha1Spa,outline=FALSE,xlab="year",ylim=c(0.95,1.2),ylab=expression(alpha[1]),main=expression(paste("Evolution of ", alpha[1], " for Spanish")))
boxplot(alpha2Spa,outline=FALSE,ylim=c(1,3),xlab="year",ylab=expression(alpha[2]),main=expression(paste("Evolution of ",alpha[2]," for Spanish")))
boxplot(xminSpa,outline=FALSE,xlab="year",ylab="xmin",main=expression("Evolution of xmin for Spanish"))
boxplot(kinkSpa,outline=FALSE,xlab="year",ylim=c(7,8),ylab="Kink",main=expression("Evolution of the kink for Spanish"))

relativeKerGer<-as.data.frame(apply(finalResultsGer,2,separatingResults,n=1))
alpha1Ger<-as.data.frame(apply(finalResultsGer,2,separatingResults,n=2))
alpha2Ger<-as.data.frame(apply(finalResultsGer,2,separatingResults,n=3))
xminGer<-as.data.frame(apply(finalResultsGer,2,separatingResults,n=4))
kinkGer<-as.data.frame(apply(finalResultsGer,2,separatingResults,n=5))
kinkGer<-as.data.frame(apply(finalResultsGer,2,separatingResults,n=5))
boxplot(relativeKerGer,outline=FALSE,xlab="year",ylab="Percentage of vocabulary in the kernel lexicon",main="Evolution of the relative size of the kernel lexicon for German",names=as.character(seq(1800,2007,1)))
boxplot(alpha1Ger,outline=FALSE,xlab="year",ylim=c(0.9,1.1),ylab="Alpha1",main="Evolution of alpha1 for German")
boxplot(alpha2Ger,outline=FALSE,ylim=c(1,3),xlab="year",ylab="Alpha2",main="Evolution of alpha2 for German")
boxplot(xminGer,outline=FALSE,xlab="year",ylab="xmin",main="Evolution of the position until which the power law holds for German")
boxplot(kinkGer,outline=FALSE,xlab="year",ylim=c(3,8),ylab="Kink",main="Evolution of the kink for German")



save.image("fixFre.RData")
