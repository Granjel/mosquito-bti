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
  cbind(alive, dead) ~ bti * food * s_time + (s_time | jar),
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
theme_set(theme_bw() + theme(strip.background = element_rect(fill = "white")))

lt <- LT50_table %>%
  mutate(
    bti = factor(bti, levels = c("No", "Yes")),
    food = factor(food, levels = c("Low", "High"))
  )

p_lt <- ggplot(lt, aes(x = LT50_hours, y = bti, color = bti)) +
  geom_errorbarh(
    aes(xmin = LT50_h_low, xmax = LT50_h_high),
    height = 0.18,
    linewidth = 0.9
  ) +
  geom_point(size = 2.8) +
  facet_wrap(
    ~food,
    nrow = 1,
    labeller = labeller(
      food = c(Low = "Low organic matter", High = "High organic matter")
    )
  ) +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

# tidy log-hours axis: auto limits from CIs, clean breaks
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
    labels = label_number(accuracy = 1),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(x = "Median survival time, LT\u2085\u2080 (hours)", y = "Bti")

# ---- 6) save -------------------------------------------------------------
ggsave(
  "figures/fig_LT50_facets.png",
  p_lt,
  width = 7.2,
  height = 3.8,
  dpi = 640
)

# ---- 7) caption suggestion -----------------------------------------------
# "Points show model-derived median survival time (LT50) per treatment;
#  whiskers are 95% CIs from a parametric bootstrap of the GLMM fixed effects.
#  Panels show low vs high organic matter; colors indicate Bti (No/

# =========================================================================
# LT50 from GLMM + figure + p-values + tables (end-to-end)
# Outputs:
#   figures/fig_LT50_facets.{pdf,png}
#   tables/LT50_by_treatment.csv
#   tables/LT50_contrasts.csv
# =========================================================================

# =========================================================================
# FULL PIPELINE: GLMM -> LT50 (CI) -> Figure (log-hours) -> Pairwise Stats
# Outputs:
#   figures/fig_LT50_facets.{pdf,png}
#   tables/LT50_by_treatment.csv
#   tables/LT50_pairwise_ratios.csv
#   tables/LT50_pairwise_diffs.csv
# =========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(MASS) # mvrnorm for parametric bootstrap
  library(scales)
})

# ------------------------ 0) Load / prepare data --------------------------
source("code/01-fetch-data.R") # must create `data`
stopifnot(exists("data"))

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
  mutate(s_time = as.numeric(scale(time_h, center = TRUE, scale = TRUE)))

# back-transform mapping
mu_time <- mean(data$time_h, na.rm = TRUE)
sd_time <- sd(data$time_h, na.rm = TRUE)

# ------------------------ 1) Fit GLMM -------------------------------------
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

# ------------------------ 2) LT50 point estimates -------------------------
# eta(s) = A + B*s; LT50 when eta=0 => s* = -A/B
mkX <- function(bti_level, food_level, s) {
  model.matrix(
    ~ bti * food * s_time,
    data = data.frame(
      bti = factor(bti_level, levels = levels(data$bti)),
      food = factor(food_level, levels = levels(data$food)),
      s_time = s
    )
  )
}

lt50_scaled_one <- function(fit, bti_level, food_level) {
  X0 <- mkX(bti_level, food_level, 0)
  X1 <- mkX(bti_level, food_level, 1)
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

# ------------------------ 3) 95% CIs via bootstrap ------------------------
set.seed(123)
Sigma <- as.matrix(vcov(m))
betas <- MASS::mvrnorm(n = 2000, mu = fixef(m), Sigma = Sigma)

X0_list <- lapply(seq_len(nrow(treat_grid)), \(i) {
  mkX(treat_grid$bti[i], treat_grid$food[i], 0)
})
X1_list <- lapply(seq_len(nrow(treat_grid)), \(i) {
  mkX(treat_grid$bti[i], treat_grid$food[i], 1)
})

lt50_boot_h <- lapply(seq_len(nrow(treat_grid)), function(i) {
  A <- as.numeric(betas %*% t(X0_list[[i]]))
  B <- as.numeric(betas %*% t(X1_list[[i]] - X0_list[[i]]))
  s_star <- -A / B
  (s_star * sd_time + mu_time)[is.finite(s_star)]
})

ci_low <- sapply(lt50_boot_h, \(x) quantile(x, 0.025, na.rm = TRUE))
ci_high <- sapply(lt50_boot_h, \(x) quantile(x, 0.975, na.rm = TRUE))

LT50_table <- treat_grid %>%
  transmute(
    bti,
    food,
    LT50_hours = round(LT50_hours, 2),
    LT50_h_low = round(ci_low, 2),
    LT50_h_high = round(ci_high, 2)
  ) %>%
  arrange(food, bti)

# ------------------------ 4) Pairwise comparisons (ALL SIX) ---------------
# aligned bootstrap draws (hours) for each treatment
build_draws <- function(bti_level, food_level) {
  X0 <- mkX(bti_level, food_level, 0)
  X1 <- mkX(bti_level, food_level, 1)
  A <- as.numeric(betas %*% t(X0))
  B <- as.numeric(betas %*% t(X1 - X0))
  s_star <- -A / B
  s_star * sd_time + mu_time
}

draws <- tibble::tibble(
  `No × Low` = build_draws("No", "Low"),
  `Yes × Low` = build_draws("Yes", "Low"),
  `No × High` = build_draws("No", "High"),
  `Yes × High` = build_draws("Yes", "High")
)

boot_ratio <- function(a, b) {
  ok <- is.finite(a) & is.finite(b) & a > 0
  r <- b[ok] / a[ok]
  c(
    est = median(r, na.rm = TRUE),
    lo = quantile(r, 0.025, na.rm = TRUE),
    hi = quantile(r, 0.975, na.rm = TRUE),
    p = 2 * pmin(mean(r >= 1, na.rm = TRUE), mean(r <= 1, na.rm = TRUE))
  )
}
boot_diff <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  d <- b[ok] - a[ok]
  c(
    est = median(d, na.rm = TRUE),
    lo = quantile(d, 0.025, na.rm = TRUE),
    hi = quantile(d, 0.975, na.rm = TRUE),
    p = 2 * pmin(mean(d >= 0, na.rm = TRUE), mean(d <= 0, na.rm = TRUE))
  )
}

labs <- colnames(draws)
pairs_idx <- combn(seq_along(labs), 2)

ratio_rows <- lapply(seq_len(ncol(pairs_idx)), function(k) {
  i <- pairs_idx[1, k]
  j <- pairs_idx[2, k]
  s <- boot_ratio(draws[[i]], draws[[j]]) # ratio = (j / i)
  tibble::tibble(
    `Comparison (B/A)` = paste(labs[j], "vs", labs[i]),
    `LT50 ratio` = as.numeric(s["est"]),
    `CI_low` = as.numeric(s["lo"]),
    `CI_high` = as.numeric(s["hi"]),
    `p` = as.numeric(s["p"])
  )
})
LT50_pairwise_ratio <- dplyr::bind_rows(ratio_rows) %>%
  mutate(
    `LT50 ratio` = scales::number(`LT50 ratio`, accuracy = 0.01),
    `95% CI` = paste0(
      "[",
      scales::number(CI_low, accuracy = 0.01),
      "–",
      scales::number(CI_high, accuracy = 0.01),
      "]"
    ),
    p = scales::pvalue(p, accuracy = 0.001)
  ) %>%
  transmute(
    `Comparison (B/A)` = `Comparison (B/A)`,
    `LT50 ratio` = `LT50 ratio`,
    `95% CI` = `95% CI`,
    `p` = p
  )

diff_rows <- lapply(seq_len(ncol(pairs_idx)), function(k) {
  i <- pairs_idx[1, k]
  j <- pairs_idx[2, k]
  s <- boot_diff(draws[[i]], draws[[j]]) # diff = (j - i)
  tibble::tibble(
    `Comparison (B−A)` = paste(labs[j], "−", labs[i]),
    `LT50 diff (h)` = as.numeric(s["est"]),
    `CI_low` = as.numeric(s["lo"]),
    `CI_high` = as.numeric(s["hi"]),
    `p` = as.numeric(s["p"])
  )
})
LT50_pairwise_diff <- dplyr::bind_rows(diff_rows) %>%
  mutate(
    `LT50 diff (h)` = scales::number(`LT50 diff (h)`, accuracy = 0.1),
    `95% CI (h)` = paste0(
      "[",
      scales::number(CI_low, accuracy = 0.1),
      "–",
      scales::number(CI_high, accuracy = 0.1),
      "]"
    ),
    p = scales::pvalue(p, accuracy = 0.001)
  ) %>%
  transmute(
    `Comparison (B−A)` = `Comparison (B−A)`,
    `LT50 diff (h)` = `LT50 diff (h)`,
    `95% CI (h)` = `95% CI (h)`,
    `p` = p
  )


# ------------------------ 5) Figure (log-hours, p-values) -----------------
# palette
if (file.exists("tools/bti-colors.R")) {
  source("tools/bti-colors.R")
} # defines `pal`
if (!exists("pal")) {
  pal <- c("No" = "#0072B2", "Yes" = "#D55E00")
} # Okabe–Ito default

theme_set(theme_bw() + theme(strip.background = element_blank()))

lt <- LT50_table %>%
  mutate(
    # display order: "No" on top; Low facet first
    bti = factor(bti, levels = c("Yes", "No")),
    food = factor(food, levels = c("Low", "High"))
  )

p_lt <- ggplot(lt, aes(x = LT50_hours, y = bti, color = bti)) +
  geom_errorbarh(
    aes(xmin = LT50_h_low, xmax = LT50_h_high),
    height = 0.18,
    linewidth = 0.9
  ) +
  geom_point(size = 2.8) +
  facet_wrap(
    ~food,
    nrow = 1,
    labeller = labeller(
      food = c(Low = "Low organic matter", High = "High organic matter")
    )
  ) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none", strip.text = element_text(face = "bold")) +
  labs(y = "Bti")

# tidy log-hours axis
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
    labels = label_number(accuracy = 1),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(x = "Median survival time, LT\u2085\u2080 (hours, log scale)")

# p-values to annotate (Yes vs No within food; and Yes High vs Low)
fmtp <- function(x) pvalue(x, accuracy = 0.001)

# compute the three p-values directly from the draws
boot_p_ratio <- function(a, b) {
  ok <- is.finite(a) & is.finite(b) & a > 0
  r <- b[ok] / a[ok]
  if (!length(r)) {
    return(NA_real_)
  }
  2 * pmin(mean(r >= 1), mean(r <= 1))
}
p_low_within <- boot_p_ratio(draws$`No × Low`, draws$`Yes × Low`)
p_high_within <- boot_p_ratio(draws$`No × High`, draws$`Yes × High`)
p_yes_between <- boot_p_ratio(draws$`Yes × Low`, draws$`Yes × High`)

ann_within <- tibble::tibble(
  food = factor(c("Low", "High"), levels = c("Low", "High")),
  x = xmax / 1.15,
  y = factor("No", levels = levels(lt$bti)), # top row
  lab = c(
    paste0("p = ", fmtp(p_low_within)),
    paste0("p = ", fmtp(p_high_within))
  )
)
ann_between_yes <- tibble::tibble(
  food = factor("High", levels = c("Low", "High")),
  x = xmax / 1.15,
  y = factor("Yes", levels = levels(lt$bti)),
  lab = paste0("Yes: High vs Low, p = ", fmtp(p_yes_between))
)

p_lt <- p_lt +
  geom_text(
    data = ann_within,
    aes(x = x, y = y, label = lab),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = -0.4,
    size = 3
  ) +
  geom_text(
    data = ann_between_yes,
    aes(x = x, y = y, label = lab),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = -0.4,
    size = 3
  )

# ------------------------ 6) Save outputs ---------------------------------
if (!dir.exists("figures")) {
  dir.create("figures", recursive = TRUE)
}
if (!dir.exists("tables")) {
  dir.create("tables", recursive = TRUE)
}

# figure
ggsave("figures/fig_LT50_facets.pdf", p_lt, width = 6.8, height = 3.6)
ggsave(
  "figures/fig_LT50_facets.png",
  p_lt,
  width = 6.8,
  height = 3.6,
  dpi = 600
)

# tables
LT50_table_pretty <- LT50_table %>%
  arrange(food, bti) %>%
  mutate(
    LT50_str = sprintf("%.2f (%.2f–%.2f)", LT50_hours, LT50_h_low, LT50_h_high)
  ) %>%
  transmute(Food = food, Bti = bti, `LT50 (hours, 95% CI)` = LT50_str)
readr::write_csv(LT50_table_pretty, "tables/LT50_by_treatment.csv")
readr::write_csv(LT50_pairwise_ratio, "tables/LT50_pairwise_ratios.csv")
readr::write_csv(LT50_pairwise_diff, "tables/LT50_pairwise_diffs.csv")

# ------------------------ Caption suggestion ------------------------------
# "Points show model-derived median survival time (LT50) per treatment;
#  whiskers are 95% CIs from a parametric bootstrap of GLMM fixed effects.
#  Panels: Low vs High organic matter; colors: Bti (No/Yes). The x-axis is
#  log-scaled (hours). Pairwise bootstrap comparisons (six total) are
#  reported in Tables Sx (ratios) and Sy (differences in hours)."
# =========================================================================
