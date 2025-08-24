# -------------------------------------------------------------------------
# model-based survival plots from a binomial GLMM (no raw summaries shown)
# outputs:
#   - fig_glmm_pred_facets.pdf   (facet by Food, color = Bti)
#   - fig_glmm_pred_onepanel.pdf (one panel, color = Bti, linetype = Food)
# -------------------------------------------------------------------------

# load packages
library(tidyverse)
library(lme4)
library(emmeans)
library(car)


# --- 0) load/prepare data -------------------------------------------------

# load data through the relevant script
source("code/01-fetch-data.R")
source("tools/bti-colors.R")

# define theme for all plots
theme_set(theme_bw() + theme(strip.background = element_blank()))

# enforce factor coding and time variables
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
    s_time = as.numeric(scale(time_h, center = TRUE, scale = TRUE))
  )

# keep the mapping to back-transform the x-axis later
mu_time <- mean(data$time_h, na.rm = TRUE)
sd_time <- sd(data$time_h, na.rm = TRUE)


# --- 1) fit or reuse the GLMM --------------------------------------------

if (!exists("m") || !inherits(m, "glmerMod")) {
  m <- glmer(
    cbind(alive, dead) ~ bti * food * s_time + (1 | experiment:unit),
    family = binomial,
    data = data,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
}

# --- diagnostics ---------------------------------------------------------

# check for overdispersion
rdf <- df.residual(m)
phi <- sum(residuals(m, type = "pearson")^2) / rdf
cat("Overdispersion ratio (Pearson/df):", round(phi, 3), "\n")

# check for singular fit
cat("isSingular:", isSingular(m, tol = 1e-4), "\n")


# --- 2) prediction grid (fixed-effects / population level) ----------------

t_seq <- seq(
  min(data$s_time, na.rm = TRUE),
  max(data$s_time, na.rm = TRUE),
  length.out = 300
)

# emmeans gives fitted means and asymptotic 95% CIs on the response scale
pred_facets <- emmip(
  m,
  bti ~ s_time | food,
  at = list(s_time = t_seq),
  type = "response",
  CIs = TRUE,
  plotit = FALSE
) %>%
  as.data.frame() %>%
  rename(pred = yvar, lower = LCL, upper = UCL) %>%
  mutate(time_h = s_time * sd_time + mu_time)


# --- 3) plots (model predictions only) ------------------------------------

# facet by Food, color = Bti
p <- ggplot(pred_facets, aes(x = time_h, y = pred, color = bti, fill = bti)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, linewidth = 0) +
  geom_line(linewidth = 1.2) +
  facet_wrap(
    ~food,
    nrow = 1,
    labeller = labeller(
      food = c(Low = "Low organic matter", High = "High organic matter")
    )
  ) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(n.breaks = 8) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Time (hours since first reading)",
    y = "Predicted survival rate",
    color = "Bti",
    fill = "Bti"
  ) +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(face = "italic")
  )


# --- 4) save files --------------------------------------------------------
ggsave(
  "figures/fig-glmm-pred-facets.jpeg",
  p,
  width = 7.2 * 0.8,
  height = 3.8 * 0.8,
  dpi = 640 * 3
)


# model diagnostics ------------------------------------------------------

# check model is fine with DHARMa
library(DHARMa)
simres <- simulateResiduals(m, n = 1000)
plot(simres)

jpeg(
  "figures/model-diagnostics/dharma-main.jpeg",
  width = 7.2 * 1.35,
  height = 3.8 * 1.35,
  res = 640 * 3,
  units = "in"
)
plot(simres)
dev.off()


# variance analysis ------------------------------------------------------

# perform Type II Wald chi-square ANOVA
anova_tab <- Anova(m, type = "II")

# convert ANOVA object to data frame
anova_df <- as.data.frame(anova_tab)

# add effect names as a proper column
anova_df$Effect <- rownames(anova_tab)

# compute significance stars based on p-values
anova_df$Signif <- symnum(
  anova_df$`Pr(>Chisq)`,
  corr = FALSE,
  na = FALSE,
  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
  symbols = c("***", "**", "*", ".", " ")
)

# format p-values: "<0.001" for very small values, else 3 decimals
anova_df$P_value <- ifelse(
  anova_df$`Pr(>Chisq)` < 0.001,
  "<0.001",
  sprintf("%.3f", anova_df$`Pr(>Chisq)`)
)

# reorder and rename columns
anova_out <- anova_df[, c("Effect", "Chisq", "Df", "P_value", "Signif")]

# export to CSV
write.csv(anova_out, "tables/anova-table.csv", row.names = FALSE)


# emmeans: survival rates + tukey ORs ------------------------------------

# define final scaled time
t_final <- max(data$s_time, na.rm = TRUE)

# emms on the link (logit) scale for contrasts (needed to get odds ratios)
emm_final_link <- emmeans(
  m,
  ~ bti * food,
  at = list(s_time = t_final)
)

# emms on the response scale for treatment-level survival probabilities
emm_final_resp <- summary(emm_final_link, type = "response")

# tukey pairwise comparisons as odds ratios (robust to column names)

# pairwise tukey contrasts on the link scale
contrasts_final <- contrast(
  emm_final_link,
  method = "pairwise",
  adjust = "tukey"
)

# back-transform to get odds ratios and 95% ci
ctr <- summary(contrasts_final, infer = c(TRUE, TRUE), type = "response")
ctr_df <- as.data.frame(ctr)

# find column names robustly across emmeans versions
or_col <- if ("odds.ratio" %in% names(ctr_df)) {
  "odds.ratio"
} else if ("ratio" %in% names(ctr_df)) {
  "ratio"
} else {
  NA_character_
}
lcl_col <- intersect(c("lower.CL", "asymp.LCL", "LCL"), names(ctr_df))[1]
ucl_col <- intersect(c("upper.CL", "asymp.UCL", "UCL"), names(ctr_df))[1]

if (is.na(or_col) || is.na(lcl_col) || is.na(ucl_col)) {
  stop(
    "could not locate OR/LCL/UCL columns. columns found: ",
    paste(names(ctr_df), collapse = ", ")
  )
}

# coerce to numeric in case they came in as character
ctr_df[[or_col]] <- as.numeric(ctr_df[[or_col]])
ctr_df[[lcl_col]] <- as.numeric(ctr_df[[lcl_col]])
ctr_df[[ucl_col]] <- as.numeric(ctr_df[[ucl_col]])

# add significance stars
ctr_df$Signif <- symnum(
  ctr_df$p.value,
  corr = FALSE,
  na = FALSE,
  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
  symbols = c("***", "**", "*", ".", " ")
)

# format p-values: "<0.001" if very small, else 3 decimals
ctr_df$p.value <- ifelse(
  ctr_df$p.value < 0.001,
  "<0.001",
  sprintf("%.3f", ctr_df$p.value)
)

# custom formatter: decimals if between 1e-3 and 1e+3, else scientific
fmt_triplet <- function(or, lcl, ucl) {
  sci <- (abs(or) < 1e-3) |
    (abs(lcl) < 1e-3) |
    (abs(ucl) < 1e-3) |
    (abs(or) >= 1e3) |
    (abs(lcl) >= 1e3) |
    (abs(ucl) >= 1e3)
  out <- character(length(or))

  # scientific notation for extreme rows
  out[sci] <- sprintf("%.2e (%.2e-%.2e)", or[sci], lcl[sci], ucl[sci])

  # decimals for the rest; use 3 decimals if any value is < 0.01 to avoid 0.00
  dec <- !sci
  need3 <- dec & ((abs(or) < 0.01) | (abs(lcl) < 0.01) | (abs(ucl) < 0.01))
  out[dec & !need3] <- sprintf(
    "%.2f (%.2f-%.2f)",
    or[dec & !need3],
    lcl[dec & !need3],
    ucl[dec & !need3]
  )
  out[need3] <- sprintf("%.3f (%.3f-%.3f)", or[need3], lcl[need3], ucl[need3])

  out
}

# apply to odds ratio and CI columns
ctr_df$OR_95CI <- fmt_triplet(
  ctr_df[[or_col]],
  ctr_df[[lcl_col]],
  ctr_df[[ucl_col]]
)

# assemble export table
contrasts_export <- ctr_df %>%
  transmute(
    Contrast = contrast,
    OddsRatio_95CI = OR_95CI,
    p_value = p.value,
    Signif = Signif
  )

# write csv
write.csv(contrasts_export, "tables/emmeans-odds-ratio.csv", row.names = FALSE)
