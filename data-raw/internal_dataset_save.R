########################################################################################
## Title: internal_data_save.R
## Author: Blanca Pijuan-Sala
## Description: Save internal data.
## Date: 30 September 2017
## **proSeq package**
########################################################################################

###SAVE INTERNAL DATASET.
#http://r-pkgs.had.co.nz/data.html
wd <- "/Users/blancap/Documents/repos_bitbucket/r-packages/00.PACKAGES/proSeq_development/proSeq/"


geneTable <- read.csv(paste0(wd,"inst/extdata/mm10_mart_export.txt"))
#save(geneTable, file="/Users/blancap/Documents/PhD_CAMBRIDGE/PhD/My_R_package/00.PACKAGES/proSeq_development/proSeq/data/geneTable.rda")

############################

setwd(paste0(wd,"R/"))
library(devtools)




devtools::use_data(geneTable, internal = TRUE,overwrite=TRUE)
