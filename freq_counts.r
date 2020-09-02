install.packages("segmented")
library(segmented)
library(ggplot2)

columns <- c("ngram", "year", "match_count", "volume_count", "language")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(getSrcDirectory()[1])

myData <- read.csv("ngrams_by_year/eng-gb-1887.csv", sep=";", col.names = columns) 
help(read.csv)


ggplot(freqCounts, aes(x=Var1, y=Freq)) + 
  geom_point() + 
  scale_y_continuous(trans='log10') + 
  scale_x_continuous(trans='log10')

freqCounts <- aggregate(data.frame(count = myData$match_count), list(frequency = myData$match_count), length)

ggplot(freqCounts, aes(x=frequency, y=count)) + 
  geom_point() + 
  scale_y_continuous(trans='log10') + 
  scale_x_continuous(trans='log10')

freqCounts$ccdf <- 1 - cumsum(freqCounts$count) / sum(freqCounts$count)

ggplot(freqCounts, aes(x=ccdf, y=frequency)) + 
  geom_point() + 
  scale_y_continuous(trans='log10') + 
  scale_x_continuous(trans='log10')

ggsave("images/eng-gb-1887-freq-counts-ccdf.png")

freqCounts <- head(freqCounts, -1)

freqCounts$logccdf <- log(freqCounts$ccdf)
freqCounts$logfreq <- log(freqCounts$frequency)


out.lm<-lm(logfreq~logccdf,data=freqCounts)

#the simplest example: the starting model includes just 1 covariate 
#.. and 1 breakpoint has to be estimated for that
o <-segmented(out.lm) #1 breakpoint for x

o$psi["psi1.logccdf", "Est."]

help(log)
exp(9.41639)



help(segmented)

o.breakpoint <- o$psi["psi1.logccdf", "Est."]

o$coefficients

o.alpha <- o$coefficients["logccdf"]
o.beta <- o$coefficients["U1.logccdf"]
o.alpha2 <- o$coefficients["U1.logccdf"] + o$coefficients["logccdf"]

o.intercept <- o$coefficients["(Intercept)"]

yy_hats <- o.intercept + o.alpha*(freqCounts$logccdf) + o.beta*pmax(freqCounts$logccdf-o.breakpoint, 0)
freqCounts["logfreq_hats"] = yy_hats

ggplot(freqCounts, aes(x=logccdf)) + 
  geom_point(aes(y=logfreq), color="blue") +
  geom_line(aes(y=logfreq_hats), color="red", size=1)

ggsave("images/breakpoint-fit-eng-gb-1887-freq-counts-ccdf.png")


