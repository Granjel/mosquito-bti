# color palette ----------------------------------------------------------

# install and/or load the package
# https://emilhvitfeldt.github.io/r-color-palettes/discrete/PNWColors/Bay/
pacman::p_load("paletteer")

# # define colors per treatment
# lvl <- c("No_Low", "Yes_Low", "No_High", "Yes_High")

# pal <- setNames(
#   paletteer_d("PNWColors::Bay")[c(2, 1, 4, 5)],
#   lvl
# )

lvl <- c("No", "Yes")

pal <- setNames(
  paletteer_d("PNWColors::Bay")[c(2, 4)],
  lvl
)
