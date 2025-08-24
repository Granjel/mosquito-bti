# =========================================================================
# LT50 from GLMM + publication-ready LT50 figure (log-hours x-axis)
# Outputs:
#   figures/fig_LT50_facets.pdf
#   figures/fig_LT50_facets.png
# =========================================================================

# ---- 0) packages ---------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(MASS) # mvrnorm for parametric bootstrap
  library(scales)
})

# ---- 1) load / prepare data ----------------------------------------------
# Your fetch script should create a tibble named `data`
source("code/01-fetch-data.R")
stopifnot(exists("data"))

# factors + time variables
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

# mapping to back-transform s_time -> hours
mu_time <- mean(data$time_h, na.rm = TRUE)
sd_time <- sd(data$time_h, na.rm = TRUE)

# ---- 2) fit GLMM (random slopes by jar; safe refit if singular) ----------
m <- glmer(
  cbind(alive, dead) ~ bti * food * s_time + (1 | experiment:unit),
  family = binomial,
  data = data,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

if (isSingular(m, tol = 1e-5)) {
  message("Note: singular fit detected; refitting with independent REs (||).")
  m <- update(m, . ~ . - (s_time | jar) + (1 + s_time || jar))
}

# ---- 3) compute LT50 from fixed effects (population level) ---------------
# eta(s) = A + B*s; LT50 when eta=0 => s* = -A/B (then back-transform to hours)
lt50_scaled_one <- function(fit, bti_level, food_level) {
  X0 <- model.matrix(
    ~ bti * food * s_time,
    data = data.frame(
      bti = factor(bti_level, levels = levels(data$bti)),
      food = factor(food_level, levels = levels(data$food)),
      s_time = 0
    )
  )
  X1 <- model.matrix(
    ~ bti * food * s_time,
    data = data.frame(
      bti = factor(bti_level, levels = levels(data$bti)),
      food = factor(food_level, levels = levels(data$food)),
      s_time = 1
    )
  )
  beta <- fixef(fit)
  A <- drop(X0 %*% beta)
  B <- drop((X1 - X0) %*% beta)
  if (!is.finite(A) || !is.finite(B) || abs(B) < .Machine$double.eps) {
    return(NA_real_)
  }
  -A / B
}

treat_grid <- expand.grid(
  bti = levels(data$bti),
  food = levels(data$food),
  KEEP.OUT.ATTRS = FALSE
)

treat_grid$LT50_s <- mapply(
  function(b, f) lt50_scaled_one(m, b, f),
  treat_grid$bti,
  treat_grid$food
)
treat_grid$LT50_hours <- treat_grid$LT50_s * sd_time + mu_time

# ---- 4) 95% CIs via parametric bootstrap of fixed effects ----------------
set.seed(123)
Sigma <- as.matrix(vcov(m))
beta_hat <- fixef(m)
n_sims <- 2000

# design rows (speed)
X0_list <- lapply(seq_len(nrow(treat_grid)), function(i) {
  model.matrix(
    ~ bti * food * s_time,
    data = data.frame(
      bti = factor(treat_grid$bti[i], levels = levels(data$bti)),
      food = factor(treat_grid$food[i], levels = levels(data$food)),
      s_time = 0
    )
  )
})
X1_list <- lapply(seq_len(nrow(treat_grid)), function(i) {
  model.matrix(
    ~ bti * food * s_time,
    data = data.frame(
      bti = factor(treat_grid$bti[i], levels = levels(data$bti)),
      food = factor(treat_grid$food[i], levels = levels(data$food)),
      s_time = 1
    )
  )
})

betas <- MASS::mvrnorm(n = n_sims, mu = beta_hat, Sigma = Sigma)

lt50_boot_h <- lapply(seq_len(nrow(treat_grid)), function(i) {
  X0 <- X0_list[[i]]
  X1 <- X1_list[[i]]
  A <- as.numeric(betas %*% t(X0))
  B <- as.numeric(betas %*% t(X1 - X0))
  s_star <- -A / B
  h_star <- s_star * sd_time + mu_time
  h_star[is.finite(h_star)]
})

ci_low <- sapply(lt50_boot_h, function(x) quantile(x, 0.025, na.rm = TRUE))
ci_high <- sapply(lt50_boot_h, function(x) quantile(x, 0.975, na.rm = TRUE))

LT50_table <- treat_grid %>%
  transmute(
    bti,
    food,
    LT50_hours = round(LT50_hours, 2),
    LT50_h_low = round(ci_low, 2),
    LT50_h_high = round(ci_high, 2)
  ) %>%
  arrange(food, bti)

# warn if any non-finite or inverted CIs
if (any(!is.finite(LT50_table$LT50_hours))) {
  warning(
    "Non-finite LT50 in some treatments; check model fit and time scaling."
  )
}
if (any(LT50_table$LT50_h_low > LT50_table$LT50_h_high, na.rm = TRUE)) {
  warning("Some LT50 CI bounds inverted; investigate bootstrap draws.")
}

print(LT50_table)

# ---- 5) figure: LT50 dot–whisker (log-hours) -----------------------------
# palette: use your file if present, otherwise a safe default
source("tools/bti-colors.R") # should define `pal`
theme_set(theme_bw() + theme(strip.background = element_rect(fill = NA)))

lt <- LT50_table %>%
  mutate(
    bti = factor(bti, levels = c("No", "Yes")),
    food = factor(food, levels = c("Low", "High"))
  )

pos <- position_dodge(width = -0.5)

p_lt <- ggplot(lt, aes(x = LT50_hours, y = food, color = bti)) +
  geom_errorbarh(
    aes(xmin = LT50_h_low, xmax = LT50_h_high),
    height = 0.18,
    linewidth = 0.9,
    position = pos
  ) +
  geom_point(size = 2.8, position = pos) +
  # ✅ Orders set within ggplot:
  scale_y_discrete(limits = c("High", "Low"), name = "Organic matter") + # Low appears on TOP
  scale_color_manual(values = pal, breaks = c("No", "Yes"), name = "Bti") + # No above Yes in legend
  labs(x = "Median survival time, LT\u2085\u2080 (hours)") +
  theme_bw() +
  theme(legend.position = "right", legend.title = element_text(face = "italic"))

# (keep your existing log10 x-axis code below)
xmin <- min(lt$LT50_h_low, na.rm = TRUE) * 0.9
xmax <- max(lt$LT50_h_high, na.rm = TRUE) * 1.1
cand <- c(2, 3, 4, 6, 8, 12, 16, 24, 36, 48, 60, 72, 96)
brks <- cand[cand >= xmin & cand <= xmax]
if (length(brks) < 3) {
  brks <- pretty(c(xmin, xmax))
}
p_lt <- p_lt +
  scale_x_log10(
    limits = c(xmin, xmax),
    breaks = brks,
    labels = scales::label_number(accuracy = 1),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.02, 0.06))
  )

p_lt <- p_lt +
  scale_x_log10(
    limits = c(xmin, xmax),
    breaks = brks,
    labels = scales::label_number(accuracy = 1),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.02, 0.06))
  )

# ---- 6) save -------------------------------------------------------------
ggsave(
  "figures/fig-LT50-facets.png",
  p_lt,
  width = 7.2 * 0.8,
  height = 3.8 * 0.8,
  dpi = 640 * 3
)


# stats! -----------------------------------------------------------------

# 1) name bootstrap vectors by treatment combination
boot_names <- with(treat_grid, paste(food, bti, sep = "_")) # e.g. "Low_No"
names(lt50_boot_h) <- boot_names
boot_named <- lt50_boot_h

# 2) function: summarize difference between two treatments (treatment2 - treatment1)
contrast_stats <- function(key2, key1) {
  diff <- boot_named[[key2]] - boot_named[[key1]]
  data.frame(
    treatment1 = key1,
    treatment2 = key2,
    diff_median_h = median(diff, na.rm = TRUE),
    ci_low_h = quantile(diff, 0.025, na.rm = TRUE),
    ci_high_h = quantile(diff, 0.975, na.rm = TRUE),
    p_value_raw = 2 *
      pmin(mean(diff > 0, na.rm = TRUE), mean(diff < 0, na.rm = TRUE))
  )
}

# 3) generate all pairwise comparisons (choose 2 of 4 = 6)
pairs <- combn(boot_names, 2, simplify = FALSE)

LT50_contrasts_all <- dplyr::bind_rows(lapply(pairs, function(p) {
  contrast_stats(p[[2]], p[[1]])
})) %>%
  dplyr::mutate(
    # round numbers
    diff_median_h = round(diff_median_h, 2),
    ci_low_h = round(ci_low_h, 2),
    ci_high_h = round(ci_high_h, 2),
    # combine into one string "median (low-high)"
    diff_h_CI = paste0(diff_median_h, " (", ci_low_h, ", ", ci_high_h, ")"),
    # format p-values
    p_value = ifelse(
      p_value_raw < 0.001,
      "<0.001",
      sprintf("%.3f", p_value_raw)
    ),
    signif = dplyr::case_when(
      p_value_raw <= 0.001 ~ "***",
      p_value_raw <= 0.01 ~ "**",
      p_value_raw <= 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  # keep only clean columns
  dplyr::select(treatment1, treatment2, diff_h_CI, p_value, signif)

print(LT50_contrasts_all)

# 4) save to CSV
write.csv(
  LT50_contrasts_all,
  "tables/table-LT50-contrasts.csv",
  row.names = FALSE
)
