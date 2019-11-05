### TITLE - PREPARE PANC8


library(Seurat)
library(SeuratData)
# InstallData("panc8")

data("panc8")

panc8 <- SplitObject(panc8, split.by = "tech")

for(proj in names(panc8)) {
  Project(panc8[[proj]]) <- proj
  panc8[[proj]]$seurat_clusters <- panc8[[proj]]$celltype
  Idents(panc8[[proj]]) <-  panc8[[proj]]$celltype
}


dir.create("./data/panc8/")
saveRDS(panc8$fluidigmc1, file = "./data/panc8/fluidigmc1.rds")
saveRDS(panc8$smartseq2, file = "./data/panc8/smartseq2.rds")
