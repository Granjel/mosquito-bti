# experiment 1 -----------------------------------------------------------

# load packages
library(tidyverse)

# load data through the relevant script
data <- readxl::read_xlsx("data/bti-raw.xlsx", sheet = 1) %>%
  mutate(experiment = 1)

# modifications in the data set
data <- data %>%
  mutate(
    # date and time need to be fixed for R
    date = as.Date(date),
    time = as.POSIXct(paste(date, format(time, "%H:%M:%S")), tz = "UTC"),
    # order bti levels
    bti = factor(
      bti,
      levels = c("C", "H"),
      labels = c("No", "Yes")
    ),
    # order food levels
    food = factor(food, levels = c("L", "H"), labels = c("Low", "High"))
  ) %>%
  # add treatment as the combination of Bti and food
  mutate(treatment = paste(bti, food, sep = "_")) %>%
  # add survival rate, number of dead, and exitus rate
  mutate(
    survival = alive / initial,
    dead = initial - alive,
    exitus = 1 - survival
  ) %>%
  # rename time because it also has the date within
  rename(datetime = time) %>%
  # rearrange the order of the columns
  relocate(
    experiment,
    date,
    datetime,
    unit,
    bti,
    food,
    treatment,
    n,
    reading,
    initial,
    alive,
    survival,
    dead,
    exitus
  )

# get color palette
source("tools/bti-colors.R")

# define theme for all plots
theme_set(theme_bw() + theme(strip.background = element_rect(fill = NA)))

data %>%
  ggplot(aes(x = reading, y = survival, color = bti, group = unit)) +
  geom_path(linewidth = 1.15) +
  geom_point(alpha = 0.33, size = 2) +
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
    # title = "Preliminary experiment",
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

# save
ggsave(
  "figures/fig-experiment1.jpeg",
  plot = last_plot(),
  width = 7.2 * 0.8,
  height = 3.8 * 0.8,
  dpi = 640 * 3
)
