install.packages("nVennR")
install.packages("factoextra")
install.packages("tidyverse")
install.packages("magrittr")
install.packages("pheatmap")
install.packages("mixtools")
install.packages("RCurl")
install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("DNAcopy")

library(devtools)
install_github("DepMap-Analytics/CoRe")
install_github("kassambara/factoextra")
install_github("francescojm/CRISPRcleanR")
