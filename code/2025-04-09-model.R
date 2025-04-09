# analysis idea ----------------------------------------------------------

# need more replicates per treatment (6 min, 10 ideally)
library(lme4)
model <- glmer(
  cbind(dead, alive) ~ bti * food * reading + (1 | jar),
  family = binomial,
  data = data
)
