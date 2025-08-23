# color palette ----------------------------------------------------------

# install and/or load the package
# https://emilhvitfeldt.github.io/r-color-palettes/discrete/PNWColors/Bay/
pacman::p_load("paletteer")

# define colors per treatment
lvl <- c("Control_Low", "High_Low", "Control_High", "High_High")

pal <- setNames(
  paletteer_d("PNWColors::Bay")[c(1, 2, 4, 5)],
  lvl
)
