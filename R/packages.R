# CRAN packages
cran_packages <- c(
  "outbreaker2",
  "ape",
  "igraph",
  "vegan",
  "distcrete",
  "tidyverse",
  "scales",
  "epitrix",
  "purrr",
  "furrr",
  "pROC",
  "ggh4x",
  "ggtext",
  "ggnewscale",
  "patchwork",
  "cowplot",
  "colorspace",
  "LaplacesDemon"
)

install.packages(cran_packages)

# GitHub packages
remotes::install_github("CyGei/simulacr")
remotes::install_github("CyGei/o2ools")
remotes::install_github("CyGei/pipetime")