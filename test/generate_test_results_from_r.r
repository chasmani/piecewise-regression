install.packages("segmented")
library(segmented)

df <- read.csv("data/test_data_simple_1_bp.csv")

fit_lm = lm(y ~ 1 + x, data = df)  # intercept-only model
fit_segmented = segmented(fit_lm, seg.Z = ~x, npsi = 1)  # Two change points along x

summary(fit_segmented)

plot(fit_segmented)
points(df)
lines.segmented(fit_segmented)
points.segmented(fit_segmented)