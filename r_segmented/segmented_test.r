install.packages("segmented")
library(segmented)
library(ggplot2)

segmented?help

help(package="segmented")  

# 1. Generate some data using a linear model

set.seed(12)
xx<-1:100

yy<-2+1.5*pmax(xx-35,0)+rnorm(100,0,2)
dati<-data.frame(x=xx,y=yy)
out.lm<-lm(y~x,data=dati)

#the simplest example: the starting model includes just 1 covariate 
#.. and 1 breakpoint has to be estimated for that
o<-segmented(out.lm) #1 breakpoint for x

ggplot(dati, aes(x=xx, y=yy)) + geom_point()

o
psi1.x
o.breakpoint <- o$psi["psi1.x", "Est."]

coef(o)
slope(o)["slope1", "Est."]
slope(o)

o$coefficients

o.alpha <- o$coefficients["x"]
o.beta <- o$coefficients["U1.x"]
o.alpha2 <- o$coefficients["U1.x"] + o$coefficients["x"]

o.intercept <- o$coefficients["(Intercept)"]

yy_hats <- o.intercept + o.alpha*(xx) + o.beta*pmax(xx-o.breakpoint, 0)
dati["yy_hats"] = yy_hats

ggplot(dati, aes(x=xx)) + 
  geom_point(aes(y=yy)) +
  geom_line(aes(y=yy_hats))

