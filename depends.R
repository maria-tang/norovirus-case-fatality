# install and load dependencies

# check pacman package installed
install.packages(setdiff("pacman", rownames(installed.packages())))

# install packages if necessary and load
pacman::p_load(
  dplyr,
  ggplot2,
  here,
  fs,
  readr,
  pammtools,
  mgcv,
  patchwork,
  gt
)
