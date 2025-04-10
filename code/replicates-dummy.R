# design and replicates --------------------------------------------------

# load packages
library(tidyverse)
library(lme4)

# settings
replicates <- 3
simple <- TRUE # FALSE for a 3x3 design

# treatments, simple or not
if (simple) {
  bti <- c("Control", "High")
  food <- c("Low", "High")
} else {
  bti <- c("Control", "Low", "High")
  food <- c("Low", "Medium", "High")
}

# readings (hours)
reading <- c(1, 2, 3, 4, 5, 6, 24)

# replicates vector
n <- 1:replicates

# create data set
data <- expand.grid(bti, food, reading, n)
data <- data %>%
  rename(bti = Var1, food = Var2, reading = Var3, n = Var4) %>%
  arrange(bti, food, n, reading) %>%
  mutate(
    alive = sample(0:8, nrow(data), replace = TRUE),
    dead = 8 - alive,
    jar = paste(bti, food, n, sep = "_")
  )

# model ------------------------------------------------------------------
#
model <- glmer(
  cbind(dead, alive) ~ (bti * food) * reading + (1 | jar),
  family = binomial,
  data = data
)

summary(model)

cat("The total number of jars is", nrow(data))
