# retrieve data from Google Drive ----------------------------------------

# load packages
library(googledrive)
library(tidyverse)
library(readxl)

# download data
drive_find(pattern = "mosquito-bti", type = "spreadsheet") %>%
  drive_download(
    type = "xlsx",
    path = "data/mosquito-bti-raw.xlsx",
    overwrite = TRUE
  )

# load data set
data <- readxl::read_xlsx("data/mosquito-bti-raw.xlsx", sheet = 1) %>%
  mutate(
    # create id tag for each jar
    jar = paste(date, mosquito, bti, food, n, sep = "_"),
    # create a column with the combination of treatments
    treatment = paste(bti, food, sep = "_"),
    # date and time need to be fixed for R
    date = as.Date(date),
    time = as.POSIXct(paste(date, format(time, "%H:%M:%S")), tz = "UTC"),
    # order bti levels # todo change if new levels are present
    bti = factor(
      bti,
      levels = c("control", "high"),
      labels = c("Control", "High")
    ),
    # order food levels # todo change if new levels are present
    food = factor(food, levels = c("low", "high"), labels = c("Low", "High")),
    # calculate percentage of dead individuals
    perc = (dead / (dead + alive)) * 100
  ) %>%
  # rename time because it also has the date within
  rename(datetime = time) %>%
  # rearrange the order of the columns
  relocate(
    type,
    date,
    datetime,
    mosquito,
    bti,
    food,
    treatment,
    n,
    jar,
    reading,
    alive,
    dead,
    perc
  )

# save as csv
write.csv(data, "data/mosquito-bti-clean.csv", row.names = FALSE)
