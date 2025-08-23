# make some preliminary plots --------------------------------------------

# load
library(tidyverse)
library(patchwork)

# load data through the relevant script
source("code/01-fetch-data.R")
source("tools/bti-colors.R")

# define theme for all plots
theme_set(theme_minimal() + theme(strip.background = element_blank()))

# kaplan–meier-style mean ± ci plots ------------------------------------

txt_sz <- 3.5

# ---- ensure treatment is an ordered factor that matches the palette
data <- data %>%
  mutate(treatment = factor(treatment, levels = lvl))

# ---- helper: mean + 95% t-ci across jars at each reading
mean_ci <- function(y) {
  y <- y[is.finite(y)]
  n <- length(y)
  m <- mean(y)
  se <- sd(y) / sqrt(max(n, 1))
  tcrit <- qt(0.975, df = max(n - 1, 1))
  data.frame(
    y = m,
    ymin = pmax(0, m - tcrit * se),
    ymax = pmin(1, m + tcrit * se)
  )
}

# ---- main plot (mean line + ribbon ci)
p <- ggplot(
  data,
  aes(x = reading, y = survival, color = treatment, group = treatment)
) +
  stat_summary(
    fun.data = mean_ci,
    geom = "ribbon",
    aes(fill = treatment),
    alpha = 0.20,
    color = NA
  ) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  scale_color_manual(values = pal, breaks = lvl, limits = lvl) +
  scale_fill_manual(values = pal, breaks = lvl, limits = lvl) +
  labs(x = "Time (hours)", y = "Survival rate") +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0),
    text = element_text(size = txt_sz * 4)
  )

# ---- build the 2×2 matrix legend (labels via annotations, not axes)
bti_lv <- c("Control", "High")
# show Food with Low on TOP, High on BOTTOM:
food_top_to_bottom <- c("Low", "High")
food_levels_pos <- rev(food_top_to_bottom) # c("High","Low") -> y=1: High, y=2: Low

legend_df <- tidyr::expand_grid(
  bti = factor(bti_lv, levels = bti_lv),
  food = factor(food_levels_pos, levels = food_levels_pos) # reversed for positioning
) %>%
  mutate(
    treatment = factor(paste(bti, food, sep = "_"), levels = lvl),
    x = as.integer(bti),
    y = as.integer(food) # y=1 -> High (bottom), y=2 -> Low (top)
  )

# offsets to push titles away from tiles
off_top_ticks <- 0.08
off_top_title <- 0.30
off_right_ticks <- 0.08
off_right_title <- 0.34

legend_mat <- ggplot(legend_df, aes(x, y, fill = treatment)) +
  geom_tile(width = 0.98, height = 0.98, color = "grey85") +
  scale_fill_manual(values = pal, limits = lvl, guide = "none") +
  # allow drawing labels outside the panel
  coord_equal(
    xlim = c(0.5, 2.5),
    ylim = c(0.5, 2.5),
    expand = FALSE,
    clip = "off"
  ) +
  theme_void() +
  # top tick labels (Bti)
  annotate(
    "text",
    x = c(1, 2),
    y = 2.5 + off_top_ticks,
    label = bti_lv,
    size = txt_sz,
    vjust = -0.25
  ) +
  # top axis title
  annotate(
    "text",
    x = 1.5,
    y = 2.5 + off_top_title,
    label = "Bti",
    fontface = "bold",
    size = txt_sz,
    vjust = -1
  ) +
  # right tick labels (Food) -> TOP: Low (y=2), BOTTOM: High (y=1)
  annotate(
    "text",
    x = 2.5 + off_right_ticks,
    y = c(2, 1),
    label = food_top_to_bottom,
    size = txt_sz,
    angle = 270,
    vjust = -0.25
  ) +
  # right axis title
  annotate(
    "text",
    x = 2.5 + off_right_title,
    y = 1.5,
    label = "Food",
    fontface = "bold",
    size = txt_sz,
    angle = 270,
    vjust = -1
  ) +
  # room on top/right so labels do not clip
  theme(plot.margin = margin(6, 25, 6, 6))

# ---- place the matrix to the right (no overlap); tweak width as needed
p_final <- p | legend_mat
p_final <- p_final + plot_layout(widths = c(1, 0.15))

# print
p_final

# save plot
ggsave(
  filename = "figures/kaplan-meier-treatments.jpeg",
  plot = p_final,
  width = 8,
  height = 5,
  dpi = 640
)
