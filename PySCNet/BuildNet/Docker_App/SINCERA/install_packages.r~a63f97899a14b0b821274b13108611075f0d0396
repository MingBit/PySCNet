packages <- c("httr", "curl", "gh", "usethis", "devtools")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biobase")

library(devtools)
devtools::install_github("xu-lab/SINCERA", force = TRUE)