# retrieve data from Google Drive ----------------------------------------

# load packages
library(googledrive)
library(tidyverse)
library(readxl)

# download data
drive_find(pattern = "sertraline_emraas", type = "spreadsheet") %>%
  drive_download(
    type = "xlsx",
    path = "sertraline/sertraline-raw.xlsx",
    overwrite = TRUE
  )

# load data and make a few adjustments
data <- readxl::read_xlsx("sertraline/sertraline-raw.xlsx", sheet = "data") %>%
  mutate(
    date = as.Date(date, format = "%d/%m/%Y"), # change date format
    dead = (5 - alive) / 5 * 100, # dead (% of the initial number)
    dose = as.factor(dose), # dose is a factor, not a dbl
    n = as.integer(n) # change replicate from dbl to int
  )

# save species name
species <- unique(data$species)

# clean df (remove unnecessary variables)
data <- data %>%
  select(date, time, dose, n, alive, deceased, dead) # we don't need columns like species or reading


# plots ------------------------------------------------------------------

# define characteristics
font.size = 20
point.size = 3
line.width = 1.5

# set theme
theme_set(
  theme_bw() +
    theme(
      text = element_text(size = font.size),
      axis.title = element_text(size = font.size * 1.15),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
)

sertraline_plot <- data %>%
  ggplot(aes(x = time, y = dead, color = dose)) +
  facet_grid(n ~ dose) +
  geom_point(size = point.size) +
  geom_path(aes(group = interaction(n, dose)), linewidth = line.width) +
  scale_color_viridis_d(begin = 0.1, end = 0.9) +
  scale_x_continuous(
    breaks = unique(data$time),
    sec.axis = sec_axis(
      ~.,
      name = "Dose of sertraline (μg/L)",
      breaks = NULL,
      labels = NULL
    )
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(
      ~.,
      name = "Replicate",
      breaks = NULL,
      labels = NULL
    )
  ) +
  labs(x = "Elapsed time (hours)", y = "Dead individuals (%)")

# save it!
ggsave(
  sertraline_plot,
  filename = "sertraline/sertraline.png",
  dpi = 640,
  width = 12,
  height = 8,
  device = "png"
)

# same but for only 48 hours
data48 <- data %>%
  filter(time <= 48)

sertraline_plot_48h <- data48 %>%
  ggplot(aes(x = time, y = dead, color = dose)) +
  facet_grid(n ~ dose) +
  geom_point(size = point.size) +
  geom_path(aes(group = interaction(n, dose)), linewidth = line.width) +
  scale_color_viridis_d(begin = 0.1, end = 0.9) +
  scale_x_continuous(
    breaks = unique(data$time),
    sec.axis = sec_axis(
      ~.,
      name = "Dose of sertraline (μg/L)",
      breaks = NULL,
      labels = NULL
    )
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(
      ~.,
      name = "Replicate",
      breaks = NULL,
      labels = NULL
    )
  ) +
  labs(x = "Elapsed time (hours)", y = "Dead individuals (%)")

# save it!
ggsave(
  sertraline_plot_48h,
  filename = "sertraline/sertraline-48h.png",
  dpi = 640,
  width = 10,
  height = 8,
  device = "png"
)

# boxplots 48 hours
sum_data <- data %>%
  filter(time %in% c(24, 48))

sertraline_boxplot <- sum_data %>%
  ggplot(aes(x = time, y = dead, group = time, color = dose, fill = dose)) +
  facet_grid(~dose) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  scale_fill_viridis_d(begin = 0.1, end = 0.9) +
  scale_color_viridis_d(begin = 0.1, end = 0.9) +
  scale_x_continuous(
    breaks = c(24, 48),
    sec.axis = sec_axis(
      ~.,
      name = "Dose of sertraline (μg/L)",
      breaks = NULL,
      labels = NULL
    )
  ) +
  labs(
    x = "Elapsed time (hours)",
    y = "Dead individuals (%)"
  ) +
  theme(panel.grid.major.x = element_blank())

ggsave(
  sertraline_boxplot,
  filename = "sertraline/sertraline-boxplot.png",
  dpi = 640,
  width = 10,
  height = 6,
  device = "png"
)
