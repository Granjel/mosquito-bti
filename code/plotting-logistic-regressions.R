library(tidyverse)

# STEP 1 — Add panel identifier to the data
data <- data %>%
  # filter(reading != 24) %>%
  mutate(panel = interaction(food, bti, sep = "_"))

# STEP 2 — Fit logistic regression per replicate
replicate_models <- data %>%
  group_by(food, bti, n, panel) %>%
  nest() %>%
  mutate(
    model = map(
      data,
      ~ suppressWarnings(glm(
        cbind(dead, alive) ~ log10(reading),
        family = binomial,
        data = .x
      ))
    ),
    pred_data = map(
      data,
      ~ tibble(
        reading = seq(min(.x$reading), max(.x$reading), length.out = 100)
      )
    ),
    pred_data = map2(
      pred_data,
      model,
      ~ mutate(.x, perc = predict(.y, newdata = .x, type = "response") * 100)
    )
  ) %>%
  unnest(pred_data)

# STEP 3 — Fit logistic regression per panel (food × bti)
panel_models <- data %>%
  group_by(food, bti, panel) %>%
  filter(n_distinct(perc) > 1) %>%
  nest() %>%
  mutate(
    model = map(
      data,
      ~ glm(cbind(dead, alive) ~ log10(reading), family = binomial, data = .x)
    ),
    pred_data = map(
      data,
      ~ tibble(
        reading = seq(min(.x$reading), max(.x$reading), length.out = 100)
      )
    ),
    pred_data = map2(
      pred_data,
      model,
      ~ mutate(.x, perc = predict(.y, newdata = .x, type = "response") * 100)
    )
  ) %>%
  unnest(pred_data)

# STEP 5 - colors for the panels
mosquito_colors <- setNames(
  paletteer_d("beyonce::X6")[c(-1, -6)],
  c("Low_High", "Low_Control", "High_Control", "High_High")
)

# STEP 4 — Plot
p1 <- ggplot() +
  # Raw data points
  geom_point(
    data = data,
    aes(x = reading, y = perc, color = panel, group = n),
    alpha = 0.33
  ) +

  # Logistic regression per replicate
  geom_line(
    data = replicate_models,
    aes(x = reading, y = perc, color = panel, group = n),
    alpha = 0.33,
    linewidth = 0.75
  ) +

  # Global panel-level logistic regression
  geom_line(
    data = panel_models,
    aes(x = reading, y = perc, group = panel, color = panel),
    linewidth = 1.25
  ) +

  facet_grid(food ~ bti) +
  scale_color_manual(values = mosquito_colors) +
  scale_x_log10(
    breaks = c(1:6, 24),
    sec.axis = sec_axis(~., name = "BTi", breaks = NULL, labels = NULL)
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~., name = "Food", breaks = NULL, labels = NULL)
  ) +
  labs(x = "Experiment time (hours)", y = "Dead individuals (%)") +
  theme(legend.position = "none")

ggsave(
  p1,
  filename = paste0("figures/", today(), "-paths.png"),
  dpi = 320,
  width = 6,
  height = 4,
  device = "png"
)
