load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Example usage
load_or_install("readxl")
load_or_install("stringr")
load_or_install("ggplot2")
load_or_install("dplyr")
load_or_install("corrplot")
load_or_install("GGally")
load_or_install("rlang")
load_or_install("ggpmisc")