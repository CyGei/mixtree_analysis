cran_packages <- c(
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
  "santoku",
  "patchwork",
  "cowplot",
  "colorspace",
  "LaplacesDemon",
  "pROC",
  "conflicted"
)
github_packages <- c("simulacr", "mixtree", "o2ools", "pipetime") #CyGei

all_packages <- c(cran_packages, github_packages)
invisible(lapply(all_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

conflicted::conflict_prefer_all(winner = "dplyr", quiet = TRUE)
rm(
  cran_packages,
  github_packages,
  all_packages
)
