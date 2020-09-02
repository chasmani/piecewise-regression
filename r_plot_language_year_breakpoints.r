install.packages("segmented")
library(segmented)
library(ggplot2)

columns <- c("ngram", "year", "match_count", "volume_count", "language")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(getSrcDirectory()[1])

