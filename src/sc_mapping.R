### TITLE : Mapping cell subpopulations between 2 single-cell samples
### AUTHOR : Javier Perales-Paton, javier.perales@bioquant.uni-heidelberg.de
### MOTIVATION : Unsupervised cell clustering on reduced dimensionality of heterogeneous single-cell transcriptomics
###         leads to k-classes that typically resemble the distinct cell types/population in the bulk tissue/sample.
###         These classes depends on within-sample heterogeneity, thus difficult to transfer to another independent sample
###         with a different magnite of heterogeneity (e.g. whole-tissue vs FACS-based sorted cells).
### METHODOLOGY : In brief, different metrics and methods are implemented to map cells:
###   - Metrics: within-sample measurement of (dis-)similarities:
###         * specScore : specicity scores by genesorteR - see ?sortGenes for more info.
###         * dge : ranking of differentially expresed genes by -log10(pval) * sign(log2FC) from wilcox test.
###   - Methods : between-sample approach to calculate cluster (dis-)similarity
###         * nes : normalized enrichment score from a pre-ranked Gene Set Enrichment Score
###         * pearson : pearson correlation
###         * spearman : spearman correlation
###         * scmap-cell : centroid projection of one sample into individuals cells from a second sample
### LICENSE : GPL-v3

### Requirements
# see ./src/create_env.sh

### Load libraries #####
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(genesorteR))
suppressPackageStartupMessages(require(ComplexHeatmap))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(scmap))


### Get parameters #####
#--- Input variables
option_list = list(
  make_option(c("--S1"), action="store", default="./data/ifnb/STIM.rds", type='character',
              help="Path to SeuratObject from Sample #1"),
  make_option(c("--S2"), action="store", default="./data/ifnb/CTRL.rds", type='character',
              help="Path to SeuratObject from Sample #2"),
  make_option(c("--METHOD"), action="store", default="pearson", type='character',
              help="Mapping method: pearson, spearman, nes, scmap-cell"),
  make_option(c("--METRIC"), action="store", default="specScore", type='character',
              help="Mapping metric: specScore, dge, centroid"),
  make_option(c("--N"), action="store", default=NA, type='numeric',
              help="Number of features (#genes) for mapping."),
  make_option(c("--RENAME"), action="store", default=NA, type='character',
              help="Rename samples."),
  make_option(c("--OUTDIR"), action="store", default="./results/ifnb", type='character',
              help="Output directory for tables and figures."),
  make_option(c("--TSK"), action="store", default=1, type='numeric',
              help="Number of cores to use for parallelization.")
)

#--- Parse parameters
opt = parse_args(OptionParser(option_list=option_list))

# Cat the input parameters
cat("[INFO] Input parameters\n", file=stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
  assign(user_input,opt[[user_input]])
}

# Fix "NA" to NA
if(!is.numeric(N)) {
  if(N=="NA") N <- NA; 
}

# Create directory if does not exist
if(!dir.exists(OUTDIR)) dir.create(OUTDIR)

# Read SeuratObject
for(Sx in ls(pattern = "S(1|2)")) {
  Sx_fl <- get(Sx)
  if(!file.exists(Sx_fl)) stop(paste0("ERROR : File '",Sx_fl,"' does not exist."))
  
  cat(paste0("[INFO] : Reading '",Sx,"' Seurat Object from '",Sx_fl,"'\n"), file = stdout())
  assign(Sx, readRDS(Sx_fl))
}

# Get the name of the samples
if(is.na(RENAME)) RENAME <- c(Project(S1),Project(S2));


### Functions by Mapping method ######

# Remove cell population with only N cell
consistent_cellpop <- function(S,min.cell=3) {
  valid_gr <- function(S, min.cells=3) {
    cnt <- table(S$seurat_clusters)
    celltype <- names(cnt)[cnt> min.cells]
    return(celltype)
  }
  
  valids <- valid_gr(S)
  if(!all(levels(Idents(S)) %in% valids)) {
    cat(paste0("[WARN] : Some cell populations were removed because n=1 : ",
               paste(setdiff(levels(Idents(S)),valids), collapse = ","),
               "\n"),file = stdout())
    S <- S[,Idents(S) %in% valid_gr(S)]
  }
  
  return(S)
}

### METRICS ####
# Calculate the specScores
calc_specScore <- function(S, genes=NULL, cores=1) {
  if(is.null(genes)) genes <- rownames(S);
  
  # Remove cell populations with less than 1 cells
  S <- consistent_cellpop(S, min.cell = 1)
  
  sg <- sortGenes(S@assays$RNA@data[genes,], Idents(S), cores = cores)
  return(as.matrix(sg$specScore))
}

# Calculate differentially expressed genes
calc_dge <- function(S, cores=1, prior.count=0.001) {
  
  # Remove cell populations with less than 3 cells
  S <- consistent_cellpop(S, min.cell = 3)
  
  # Calculate differentially expressed genes between clusters
  cat(paste0("[INFO] Be patient...\n"),file=stdout())
  if(cores>1) {
    cl <- makeCluster(cores)
    clusterExport(cl, "S")
    dge <-  parLapply(cl, levels(S), function(idx) {
      dge.idx <- Seurat::FindMarkers(S, ident.1 = idx, ident.2 = setdiff(levels(S), idx),
                                     logfc.threshold = 0, min.pct = 0, min.diff.pct = -Inf)
    })
    stopCluster(cl)
  } else {
    dge <-  sapply(levels(S), function(idx) {
      dge.idx <- Seurat::FindMarkers(S, ident.1 = idx, ident.2 = setdiff(levels(S), idx),test.use="wilcox",
                                     logfc.threshold = 0, min.pct = 0, min.diff.pct = -Inf)
    },simplify = FALSE,USE.NAMES = TRUE)
  }
  # Get genes
  genes <- unique(unlist(lapply(dge, rownames)))
  
  # Create matrix of dge
  rnk <- lapply(dge, function(z) setNames(-log10(z[genes,"p_val"]+prior.count) * sign(z[genes,"avg_logFC"]), genes))
  rnk <- do.call("cbind",rnk)
  return(rnk)
}

### METHODS #####
# Calculate NES
calc_nes <- function(rnk,gs.list) {
  nes <- apply(rnk, 2, function(z) {
    res <- fgsea::fgsea(gs.list, stats = z, nperm = 1000)
    nes <- setNames(res$NES,res$pathway)
  })
  return(nes)
}

# Within-cluster proportion of cells projected based on centroids
# aka. scmapCell
calc_projcentr <- function() {
  
}

# Estimate distance
sc_dist <- function(S1, S2, method=c("pearson"), metric="specScore", Nfeatures=25, cores=1) {
  if (metric=="specScore") {
    common <- intersect(rownames(S1),rownames(S2))
    M1 <- calc_specScore(S1, genes = common, cores=cores)
    M2 <- calc_specScore(S2, genes = common, cores=cores)
  } else if (metric == "dge") {
    M1 <- calc_dge(S1)
    M2 <- calc_dge(S2)
  }
  
  genes <- NULL
  if(!is.na(Nfeatures)) {
    feat1 <- rownames(M1)[rowSums(apply(-M1,2,rank) <= Nfeatures) > 0]
    feat2 <- rownames(M2)[rowSums(apply(-M2,2,rank) <= Nfeatures) > 0]
    genes <- unique(c(intersect(feat1,rownames(M2)), intersect(feat2, rownames(M1))))
  }
  if(is.null(genes)) genes <- intersect(rownames(S1),rownames(S2))
  
  if (method %in% c("pearson","spearman")) {
    res <- cor(M1[genes,],M2[genes,], method = method)
  } else if (method=="nes") {
    GS1 <- sapply(colnames(M1), function(z) names(sort(M1[,z],decreasing = TRUE))[1:Nfeatures], simplify = FALSE)
    GS2 <- sapply(colnames(M2), function(z) names(sort(M2[,z],decreasing = TRUE))[1:Nfeatures], simplify = FALSE)
    
    
    NES1 <- calc_nes(M1, GS2)
    NES1 <- t(NES1)
    
    NES2 <- calc_nes(M2, GS1)
    
    res <- (NES1 + NES2)/2
    # res <- matrix(NA, nrow=nrow(NES1), ncol=ncol(NES2), dimnames=list(rownames(NES1), colnames(NES2)))
    # res[lower.tri(res, diag=FALSE)] <- NES1[lower.tri(NES1,diag = FALSE)]
    # res[upper.tri(res, diag=FALSE)] <- NES2[upper.tri(NES2,diag=FALSE)]
  } else if (method=="scmap-cell") {
    SC1 <- seurat2scmap(S1)
    SC2 <- seurat2scmap(S2)
    # Create a common space
    common <- intersect(rownames(SC1),rownames(SC2))
    # Subset to same space
    SC1 <- SC1[common,]
    SC2 <- SC2[common,]
    rm(S1,S2) # Alliviate memory
    
    SC2 <- selectFeatures(SC2, 500, TRUE)
    SC2 <- indexCell(SC2)
    
    scmapCell_results <- scmapCell(SC1, list(S1 = metadata(SC2)$scmap_cell_index))
    scmapCell_clusters <- scmapCell2Cluster(scmapCell_results, list(as.character(colData(SC2)$cell_type1)))
    # Calculate the proportion of cells that were classified as S2 cell type
    res <- table(SC1$seurat_clusters,scmapCell_clusters$scmap_cluster_labs)*100 / ncol(SC1)
    attributes(res)$class <- "matrix"
  }
  
  return(res)
}


# Dist matrix visualization
dist_vis <- function(mat, Stat="Stat", tls=c("S1","S2"), wh="heatmap") {
  if(wh=="heatmap") {
    Heatmap(mat, name=Stat,
            row_title = tls[1],row_title_gp = gpar(fontsize=32),
            column_title = tls[2], column_title_gp = gpar(fontsize=32),
            row_names_side = "left", row_names_gp = gpar(fontsize=26),
            column_names_side = "top", column_names_gp = gpar(fontsize=26),
            heatmap_legend_param= list(legend_height = unit(4, "cm"),
                                       title_gp=gpar(fontsize=16),
                                       labels_gp = gpar(fontsize = 15)),
            column_names_rot = 90, cluster_columns = FALSE, cluster_rows = FALSE) 
  } else if(wh=="graph") {
    
  }
  
}

# Pkg compatibility
seurat2scmap <- function(S) {
  SC <- as.SingleCellExperiment(S)
  SC$cell_type1 <- S$seurat_clusters
  rowData(SC)$feature_symbol <- rownames(SC)
  # This is very bad for MEM performance... but scmap assumes normal matrices but not sparse
  logcounts(SC) <- as.matrix(logcounts(SC))
  counts(SC) <- as.matrix(counts(SC))
  
  return(SC)
}

### Main #####
res <- sc_dist(S1,S2,method=METHOD, metric=METRIC,Nfeatures = N, cores=TSK)

## Save simil matrix
write.table(res,file = paste0(OUTDIR,"/","res_",paste(RENAME,collapse = "-"),"_",METHOD,"_",METRIC,"_",N,".tsv"),
            sep="\t", quote = FALSE, row.names = TRUE, col.names = NA)

## Visual plot
# Heatmap
png(paste0(OUTDIR,"/","heatmap_",paste(RENAME,collapse = "-"),"_",METHOD,"_",METRIC,"_",N,".png"),
    height = 800*3, width = 800*3, res=280)
print(dist_vis(res, METHOD, tls = RENAME, wh="heatmap"))
dev.off()

# Graph
#png(paste0(OUTDIR,"/","graph_",paste(RENAME,collapse = "-"),"_",METHOD,"_",METRIC,"_",N,".png"),
#    height = 800*3, width = 800*3, res=280)
#print(dist_vis(res, METHOD, tls = RENAME, wh="heatmap"))
#dev.off()

### Log run ####
sink(file = paste0(OUTDIR,"/","sessionInfo_",paste(RENAME,collapse = "-"),"_",METHOD,"_",METRIC,"_",N,".txt"))
## Parameters used
cat("[INFO] Input parameters\n", file=stdout())
for(user_input in names(opt)) {
  if(user_input=="help") next;
  cat(paste0("[INFO] ",user_input," => ",opt[[user_input]],"\n"),file = stdout())
}
## Pkg version
sessionInfo()
sink()

