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




wd <- "/Users/blancap/Documents/PhD_CAMBRIDGE/PhD/Experiments/"

counts = read.table(paste0(wd, "PhD_BPS13/RUN2_B6_E825/data/PhD_BPS13_run2_countsRaw.txt"),sep="\t",header=TRUE)

colnames(counts)<-sub("SLX-","SLX.",gsub(".H.*","",colnames(counts)))
colnames(counts)<-sub("_",".",colnames(counts))
#library(xlsx)

#metadata <- read.xlsx(paste0(wd,"PhD_BPS13/RUN2_B6_E825/data/PhD_BPS13_run2_metadata.xlsx"),sheetIndex = 1)


#metadata <- metadata[,is.na((metadata[1,]))==FALSE]
#metadata <- metadata[is.na((metadata[,1]))==FALSE,]

#save(metadata,file=paste0(wd,"PhD_BPS13/RUN2_B6_E825/data/PhD_BPS13_run2_metadata.rda"))
load(paste0(wd,"PhD_BPS13/RUN2_B6_E825/data/PhD_BPS13_run2_metadata.rda"))
rownames(metadata) <- as.character(paste0(metadata$CRI.identifier,".",metadata$CI.annotation.of.index))
rownames(metadata)<-sub("-",".",colnames(counts))


metadata1 <- metadata[,c("Plate.number","Position.in.96.well.plate_sorted","Cell.type...general")]
colnames(metadata1)<-c("PlateNumber","well","celltype")

counts0 <- counts
countsMatrix <- counts0[,(metadata1$celltype)%in%c("YolkSac","Allantois")]



meta <- metadata1[(metadata1$celltype)%in%c("YolkSac","Allantois"),]
meta$celltype <- as.character(meta$celltype)
meta$celltype[meta$celltype=="YolkSac"] <- "YS"
meta$celltype[meta$celltype=="Allantois"] <- "AL"
s <- sample(rownames(meta),40)
meta <- meta[s,]
countsMatrix <- countsMatrix[,s]
wd <- "/Users/blancap/Documents/repos_bitbucket/r-packages/00.PACKAGES/proSeq_development/proSeq/"
save(meta, file=paste0(wd,"data/meta.rda"))
save(countsMatrix, file=paste0(wd,"data/countsMatrix.rda"))
