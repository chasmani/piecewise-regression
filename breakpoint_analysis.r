install.packages("segmented")
library(segmented)
library(ggplot2)


columns <- c("ngram", "year", "match_count", "volume_count", "language")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(getSrcDirectory()[1])

myData <- read.csv("ngrams_by_year/eng-gb-1887.csv", sep=";", col.names = columns) 
help(read.csv)

# Get over 20 matches only
myData <- subset(myData, match_count > 20)


# Set ranks
myData <- myData[order(myData$match_count, decreasing=TRUE),]
myData$rank = 1:nrow(myData)

myData$ccdf <- 1 - cumsum(myData$match_count) / sum(myData$match_count)
myData$logccdf <- log(myData$ccdf)



#
ggplot(myData, aes(x=rank, y=ccdf)) + 
  geom_point() + 
  scale_y_continuous(trans='log10') + 
  scale_x_continuous(trans='log10')

ggsave("images/eng-gb-1887-pdf.png")


myData$logrank <- log(myData$rank)
myData$logfreq <- log(myData$match_count)

ggplot(myData, aes(x=logrank, y=logfreq)) + 
  geom_point()


out.lm<-lm(logfreq~logrank,data=myData)

#the simplest example: the starting model includes just 1 covariate 
#.. and 1 breakpoint has to be estimated for that
o <-segmented(out.lm, npsi=1) #1 breakpoint for x
help(segmented)

o$psi["psi1.logrank", "Est."]

help(log)
exp(9.41639)



help(segmented)

o.breakpoint <- o$psi["psi1.logrank", "Est."]

o$coefficients

o.alpha <- o$coefficients["logrank"]
o.beta <- o$coefficients["U1.logrank"]
o.alpha2 <- o$coefficients["U1.logrank"] + o$coefficients["logrank"]

o.intercept <- o$coefficients["(Intercept)"]

yy_hats <- o.intercept + o.alpha*(myData$logrank) + o.beta*pmax(myData$logrank-o.breakpoint, 0)
myData["logfreq_hats"] = yy_hats

ggplot(myData, aes(x=logrank)) + 
  geom_point(aes(y=logfreq), color="blue") +
  geom_line(aes(y=logfreq_hats), color="red", size=1)

ggplot(myData, aes(x=logrank)) + 
  geom_point(aes(y=logfreq))

ggplot(myData, aes(x=logrank)) + 
  geom_line(aes(y=logfreq_hats))

ggsave("images/breakpoint_fit_eng-gb-1887.png")



## Slightly more complicated - 2 breakpoints
o2 <-segmented(out.lm, npsi=2, psi=c(9, 10)) #1 breakpoint for x

o$psi

o2$psi
o2$psi["psi1.logrank", "Est."]

o2.breakpoint1 <- o2$psi["psi1.logrank", "Est."]
o2.breakpoint2 <- o2$psi["psi2.logrank", "Est."]

o2$coefficients
o2.alpha <- o2$coefficients["logrank"]
o2.beta1 <- o2$coefficients["U1.logrank"]
o2.beta2 <- o2$coefficients["U2.logrank"]
o2.intercept <- o2$coefficients["(Intercept)"]


yy_hats <- o2.intercept + o2.alpha*(myData$logrank) + o2.beta1*pmax(myData$logrank-o2.breakpoint1, 0) + o2.beta2*pmax(myData$logrank-o2.breakpoint2, 0)

myData["logfreq_hats"] = yy_hats

ggplot(myData, aes(x=logrank)) + 
  geom_point(aes(y=logfreq), color="green", size=1) +
  geom_line(aes(y=logfreq_hats), color="purple", size=1, linetype = "dotted") +
  geom_vline(xintercept=o2.breakpoint1, color="purple") + 
  geom_vline(xintercept=o2.breakpoint2, color="purple")

ggsave("images/breakpoint-double-eng-gb-1887-pdf.png")


exp(8.69)
exp(10.1)