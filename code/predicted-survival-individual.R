# =========================================================================
# GLMM predicted survival (population curves + jar-specific predicted lines)
# Output: figures/fig-glmm-pred-facets-individual.jpeg
# =========================================================================

# ---- packages ------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(emmeans)
})

# ---- 0) load / prepare data ---------------------------------------------
source("code/01-fetch-data.R") # must create `data`
source("tools/bti-colors.R")

# theme
theme_set(theme_bw() + theme(strip.background = element_rect(fill = NA)))

# factors + time vars
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

# mapping to back-transform (for labels if needed)
mu_time <- mean(data$time_h, na.rm = TRUE)
sd_time <- sd(data$time_h, na.rm = TRUE)

# ---- 1) fit or reuse GLMM -----------------------------------------------
if (!exists("m") || !inherits(m, "glmerMod")) {
  m <- glmer(
    cbind(alive, dead) ~ bti * food * s_time + (s_time | jar),
    family = binomial,
    data = data,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
  if (isSingular(m, tol = 1e-5)) {
    message("Note: singular fit; refitting with independent REs (||).")
    m <- update(m, . ~ . - (s_time | jar) + (1 + s_time || jar))
  }
}

# ---- 2) population-level predictions (fixed effects only) ----------------
t_seq <- seq(
  min(data$s_time, na.rm = TRUE),
  max(data$s_time, na.rm = TRUE),
  length.out = 300
)

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
  mutate(
    time_h = s_time * sd_time + mu_time,
    food = factor(food, levels = c("Low", "High")),
    bti = factor(bti, levels = c("No", "Yes"))
  )

# ---- 2b) jar-level predicted trajectories (include random effects) -------
data_pred_jar <- data %>%
  mutate(
    pred_jar = predict(m, newdata = ., type = "response", re.form = NULL)
  )

# ---- 3) plot: jar predicted lines under population ribbons/curves --------
p <- ggplot() +
  # (i) thin jar-specific predicted lines
  geom_line(
    data = data_pred_jar,
    aes(x = time_h, y = pred_jar, group = jar, color = bti),
    linewidth = 0.5,
    alpha = 0.25
  ) +
  # (ii) population-level 95% CI ribbons
  geom_ribbon(
    data = pred_facets,
    aes(x = time_h, ymin = lower, ymax = upper, fill = bti),
    alpha = 0.18,
    linewidth = 0
  ) +
  # (iii) raw data points
  geom_point(
    data = data,
    aes(x = time_h, y = survival, color = bti),
    alpha = 0.2,
    size = 1
  ) +
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
    y = "Survival rate",
    color = "Bti",
    fill = "Bti"
  ) +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(face = "italic")
  )

# ---- 4) save -------------------------------------------------------------
ggsave(
  "figures/fig-glmm-pred-facets-individual.jpeg",
  p,
  width = 7.2 * 0.8,
  height = 3.8 * 0.8,
  dpi = 640 * 3
)

# ---- caption snippet -----------------------------------------------------
# "Thin lines show jar-specific predicted survival (GLMM, including random effects).
#  Thick lines and ribbons show population-level predictions (GLMM fixed effects)
#  with asymptotic 95% CIs from emmeans. Panels: Low vs High organic matter;
#  colors: Bti (No/Yes)."
# =========================================================================
