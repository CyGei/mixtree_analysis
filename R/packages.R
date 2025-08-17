if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(
  outbreaker2,
  ape,
  igraph,
  vegan,
  distcrete,
  tidyverse,
  scales,
  epitrix,
  purrr,
  furrr,
  pROC,
  ggh4x,
  ggtext,
  ggnewscale,
  patchwork,
  cowplot,
  colorspace,
  LaplacesDemon
)
pacman::p_load_gh("CyGei/simulacr") # to simulate outbreaks
pacman::p_load_gh("CyGei/o2ools") # helper functions for outbreaker2
pacman::p_load_gh("CyGei/pipetime") # time |> operations
