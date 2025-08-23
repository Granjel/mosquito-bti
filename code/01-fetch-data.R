# retrieve data from Google Drive ----------------------------------------

# load packages
library(googledrive)
library(tidyverse)
library(readxl)

# download data
drive_find(pattern = "bti-data", type = "spreadsheet") %>%
  drive_download(
    type = "xlsx",
    path = "data/final-bti-raw.xlsx",
    overwrite = TRUE
  )

# load data set from different experiments and label them
exp2 <- readxl::read_xlsx("data/final-bti-raw.xlsx", sheet = 2) %>%
  mutate(experiment = 2)
exp3 <- readxl::read_xlsx("data/final-bti-raw.xlsx", sheet = 3) %>%
  mutate(experiment = 3)

# merge them in a unique data set
data <- rbind(exp2, exp3) %>% relocate(experiment)

# modifications in the data set
data <- data %>%
  mutate(
    # create id tag for each jar
    jar = paste(experiment, unit, bti, food, n, sep = "_"),
    # date and time need to be fixed for R
    date = as.Date(date),
    time = as.POSIXct(paste(date, format(time, "%H:%M:%S")), tz = "UTC"),
    # order bti levels
    bti = factor(
      bti,
      levels = c("C", "H"),
      labels = c("Control", "High")
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
    jar,
    reading,
    initial,
    alive,
    survival,
    dead,
    exitus
  )

# clean environment
rm(exp2, exp3)

# save as csv
write.csv(data, "data/bti-clean.csv", row.names = FALSE)
