cran_packages <- c(
  "outbreaker2",
  "ape",
  "igraph",
  "vegan",
  "distcrete",
  "tidyverse",
  "arrow",
  "scales",
  "epitrix",
  "purrr",
  "furrr",
  "pROC",
  "ggh4x",
  "ggtext",
  "ggnewscale",
  "tidygraph",
  "ggraph",
  "patchwork",
  "cowplot",
  "colorspace",
  "LaplacesDemon"
)
github_packages <- c("simulacr", "o2ools", "pipetime") #CyGei

all_packages <- c(cran_packages, github_packages)
invisible(lapply(all_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

rm(
  cran_packages,
  github_packages,
  all_packages
)
