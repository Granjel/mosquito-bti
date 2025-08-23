# analysis idea ----------------------------------------------------------

library(lme4)
model <- glmer(
  cbind(dead, alive) ~ bti * food * reading + (1 | experiment) + (1 | unit),
  family = binomial,
  data = data
)


summary(model)
