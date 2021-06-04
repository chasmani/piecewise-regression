install.packages("segmented")
library(segmented)
library(ggplot2)

# PLant has breakpoints but in seperate groups
data(plant)
ggplot(plant, aes(x=time, y=y)) + geom_point(shape=plant$group)

# Down doesn't seem to have a breakpoint
data(down)
ggplot(down, aes(x=age, y=cases)) + geom_point()

# Stagnant has some basic breakpoint data
data(stagnant)
ggplot(stagnant, aes(x=x, y=y)) + geom_point()

out.lm <-lm(y~x, data=stagnant)
o <- segmented(out.lm)

draw.history(o, "x")
plot.segmented(o)
points(o)

# Simulate
set.seed(42)  # I always use 42; no fiddling
df = data.frame(
  x = 1:100,
  y = c(rnorm(30, 2), rnorm(40, 0), rnorm(30, 1))
)

# Plot it
plot(df)

fit_lm = lm(y ~ 1 + x, data = df)  # intercept-only model
fit_segmented = segmented(fit_lm, seg.Z = ~x, npsi = 2)  # Two change points along x

summary(fit_segmented)

plot(fit_segmented)
points(df)
lines.segmented(fit_segmented)
points.segmented(fit_segmented)








