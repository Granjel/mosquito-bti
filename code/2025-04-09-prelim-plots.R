# make some preliminary plots --------------------------------------------

# load data through the relevant script
source("code/01-fetch-data.R")

# load colors
source("tools/mosquito-colors.R")

# define theme for all plots
theme_set(theme_bw() + theme(strip.background = element_blank()))


# paths in time ----------------------------------------------------------

p1 <- data %>%
  ggplot(aes(x = reading, y = perc, color = bti, group = n)) +
  facet_grid(
    food ~ bti
  ) +
  geom_point() +
  geom_path() +
  labs(x = "Experiment time (hours)", y = "Dead individuals (%)") +
  scale_color_manual(values = mosquito_colors) +
  scale_y_continuous(
    sec.axis = sec_axis(~., name = "Food", breaks = NULL, labels = NULL)
  ) +
  scale_x_log10(
    breaks = c(1:6, 24),
    sec.axis = sec_axis(~., name = "BTi", breaks = NULL, labels = NULL)
  ) +
  theme(legend.position = "none")

ggsave(
  p1,
  filename = paste0("figures/", today(), "-paths.png"),
  dpi = 320,
  width = 6,
  height = 4,
  device = "png"
)

# barplots ---------------------------------------------------------------

p2 <- data %>%
  ggplot(aes(x = food, y = perc, fill = bti)) +
  facet_grid(
    . ~ bti,
    labeller = labeller(
      # todo change if more levels are added
      bti = c("Control" = "BTi: Control", "High" = "BTi: High")
    )
  ) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Food", y = "Dead individuals (%)") +
  scale_fill_manual(values = mosquito_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(legend.position = "none")

ggsave(
  p2,
  filename = paste0("figures/", today(), "-barplot.png"),
  dpi = 320,
  width = 5,
  height = 3,
  device = "png"
)
