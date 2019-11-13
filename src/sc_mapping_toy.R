### TITLE : [TOY script] Mapping cell subpopulations between 2 single-cell samples
### AUTHOR : Javier Perales-Paton, javier.perales@bioquant.uni-heidelberg.de
### LICENSE : GPL-v3
### MANIFEST #### 
# Minimal script of'https://github.com/saezlab/sc_map/blob/master/src/sc_mapping.R',
#   incl. only one module for one of the methods available.
# Input data: run 'https://github.com/saezlab/sc_map/blob/master/src/prepare_panc8.R'

### Load libraries #####
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(ComplexHeatmap))

### Input parameters #####
S1 <- "./data/panc8/fluidigmc1.rds"
S2 <- "./data/panc8/smartseq2.rds"
N <- 2000
OUTDIR <- "./results/panc8"
# Create directory if does not exist
if(!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# Read the two SeuratObjects
for(Sx in ls(pattern = "S(1|2)")) {
  Sx_fl <- get(Sx)
  if(!file.exists(Sx_fl)) stop(paste0("ERROR:File '",Sx_fl,"' does not exist."))
  cat(paste0("Reading '",Sx,"' Object from '",Sx_fl,"'\n"), file = stdout())
  assign(Sx, readRDS(Sx_fl))
}

### FUNCTIONS ######
# Select features
sel_features <- function(S, N=2000) {
  S <- NormalizeData(object = S, verbose = FALSE)
  S <- FindVariableFeatures(object = S, selection.method = "vst", nfeatures = N)
  sel <- VariableFeatures(S)
  
  return(sel)
}

get_anchored_preds <-function(Squery, Sref, Sref.features) {
  Squery.anchors <- FindTransferAnchors(reference = Sref, query=Squery,
                                        features = Sref.features, dims=1:30)
  Squery.preds <- TransferData(anchorset = Squery.anchors,
                               refdata = Sref$seurat_clusters, dims = 1:30)
  Squery.predids <- setNames(Squery.preds$predicted.id, rownames(Squery.preds)) 
  
  stopifnot(all(colnames(Squery) == names(Squery.predids)))
  return(Squery.predids)
}

### METRICS ####
calc_prop <- function(S, cell_class) {
  S_cmat <- table(S$seurat_clusters,cell_class)
  S_totalclass <- table(S1$seurat_clusters)
  res <- sweep(S_cmat, 1, S_totalclass, "/")
  
  if(is(res)[1]=="table") attributes(res)$class <- "matrix";
  
  return(res)
}

### METHODS #####
# Estimate distance
sc_dist <- function(S1, S2, Nfeatures=25, cores=1) {
  # Select number of features
  genes <- NULL
  if(!is.na(Nfeatures)) { # N seleced genes
    genes <- sel_features(S2, N=Nfeatures)
  } else { # All genes
    genes <- intersect(rownames(S1),rownames(S2))
  }
  
  # Transfer labels
  S1_preds <- get_anchored_preds(Squery = S1, Sref = S2, Sref.features = genes)
  res <- calc_prop(S = S1, cell_class = S1_preds)
  
  return(res)
}

# Dist matrix visualization
dist_vis <- function(mat, Stat="Prop cells", tls=c("S1","S2")) {
    color_fun <- circlize::colorRamp2(c(0,+1), c("white","red"))
    Heatmap(mat, name=Stat,
            col=color_fun,
            row_title = tls[1],row_title_gp = gpar(fontsize=32),
            column_title = tls[2], column_title_gp = gpar(fontsize=32),
            row_names_side = "left", row_names_gp = gpar(fontsize=26),
            column_names_side = "top", column_names_gp = gpar(fontsize=26),
            heatmap_legend_param= list(legend_height = unit(4, "cm"),
                                       title_gp=gpar(fontsize=16),
                                       labels_gp = gpar(fontsize = 15)),
          column_names_rot = 90, cluster_columns = FALSE, cluster_rows = FALSE) 
}

### MAIN #####
# Get simil matrix (proportion) between clusters from two independent samples
res <- sc_dist(S1,S2,Nfeatures = N)
## Visualize simil matrix
draw(dist_vis(res, tls = c(Project(S1), Project(S2))),
     column_title=paste0("Nfeatures=",N), column_title_gp=gpar(fontsize=44))
