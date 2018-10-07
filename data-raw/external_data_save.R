########################################################################################
## Title: external_data_save.R
## Author: Blanca Pijuan-Sala
## Description: Save external data.
## Date: 30 September 2017
## **proSeq package**
########################################################################################

###SAVE EXTERNAL DATASET.
#http://r-pkgs.had.co.nz/data.html

wd <- "/Users/blancap/Documents/PhD_CAMBRIDGE/PhD/My_R_package/00.PACKAGES/proSeq_development/proSeq/"


geneTable <- read.csv(paste0(wd,"inst/extdata/mm10_mart_export.txt"))
save(geneTable, file=paste0(wd,"data/geneTable.rda"))
