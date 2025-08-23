# differences between treatments -----------------------------------------

# load
library(tidyverse)
library(lme4)
library(emmeans)
library(car)

# load data through the relevant script
source("code/01-fetch-data.R")

# add random variable
data <- data %>% mutate(random = paste(experiment, unit, sep = "_"))


# wow --------------------------------------------------------------------

# load packages
library(lme4)
library(car)
library(emmeans)

# prep data (assuming you've already run your 01-fetch-data.R)
data <- data %>%
  mutate(
    bti = factor(bti, levels = c("No", "Yes")),
    food = factor(food, levels = c("Low", "High"))
  ) %>%
  group_by(jar) %>%
  mutate(
    time_h = as.numeric(difftime(datetime, min(datetime), units = "hours"))
  ) %>%
  ungroup() %>%
  mutate(
    s_time = as.numeric(scale(time_h, center = TRUE, scale = TRUE)),
    obs = row_number()
  )

# fit glmm with time × treatments, random slopes per jar
m <- glmer(
  cbind(alive, dead) ~ bti * food * s_time + (s_time || experiment / unit),
  family = binomial,
  data = data,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

summary(m)

# check overdispersion
overdisp_fun <- function(model) {
  rp <- residuals(model, type = "pearson")
  n <- length(rp)
  p <- length(fixef(model)) + length(VarCorr(model))
  rdf <- n - p
  pearson.chisq <- sum(rp^2)
  ratio <- pearson.chisq / rdf
  pval <- pchisq(pearson.chisq, df = rdf, lower.tail = FALSE)
  data.frame(chisq = pearson.chisq, ratio = ratio, rdf = rdf, p = pval)
}
overdisp_fun(m)

# if ratio > ~1.2, refit with obs-level random effect
if (overdisp_fun(m)$ratio[1] > 1.2) {
  m <- update(m, . ~ . + (1 | obs))
}

# test fixed effects (wald chi-square)
Anova(m, type = "II")

# treatment contrasts at the final time point
t_final <- max(data$s_time)
emm_final <- emmeans(
  m,
  ~ bti * food,
  at = list(s_time = t_final),
  type = "response"
)
pairs(emm_final, adjust = "tukey")
