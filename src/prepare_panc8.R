### TITLE : Prepare 'PANC8' for sc_mapping as a toy dataset
### AUTHOR: Javier Perales-Paton, javier.perales@bioquant.uni-heidelberg.de
### LICENSE: GPL-v3

library(Seurat)
library(SeuratData)
# InstallData("panc8") # Just the 1st time, wrapper to download the data

data("panc8")

panc8 <- SplitObject(panc8, split.by = "tech")

for(proj in names(panc8)) {
  Project(panc8[[proj]]) <- proj
  panc8[[proj]]$seurat_clusters <- panc8[[proj]]$celltype
  Idents(panc8[[proj]]) <-  panc8[[proj]]$celltype
}

if(!dir.exists("./data/panc8")) dir.create("./data/panc8/", recursive = TRUE)
saveRDS(panc8$fluidigmc1, file = "./data/panc8/fluidigmc1.rds")
saveRDS(panc8$smartseq2, file = "./data/panc8/smartseq2.rds")