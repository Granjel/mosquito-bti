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
    cbind(alive, dead) ~ bti * food * s_time + (s_time || experiment / unit),
    family = binomial,
    data = data,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
}


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
    y = "Predicted survival rate (GLMM)",
    color = "Bti",
    fill = "Bti"
  ) +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  )


# --- 4) save files --------------------------------------------------------
ggsave(
  "figures/fig_glmm_pred_facets.jpeg",
  p,
  width = 7.2,
  height = 3.8,
  dpi = 640
)

# --- 5) suggested caption snippet ----------------------------------------
# "Lines show predicted survival from the binomial GLMM (fixed-effects / population level).
#  Shaded ribbons are 95% confidence intervals from the model (emmeans)."
