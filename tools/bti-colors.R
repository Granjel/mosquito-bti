# color palette ----------------------------------------------------------

# install and/or load the package
# https://emilhvitfeldt.github.io/r-color-palettes/discrete/PNWColors/Bay/
pacman::p_load("paletteer")

# define colors per treatment
bti_colors <- setNames(
  paletteer_d("PNWColors::Bay")[c(2, 4)],
  c("Control", "High")
) # todo change if more levels are added
