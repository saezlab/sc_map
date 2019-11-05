```sh
conda create -n scmap r-base=3.6.1 bioconductor-complexheatmap r-seurat parallel r-optparse
conda install -n scmap bioconductor-scmap
conda install -n scmap tar zip

conda activate scmap
$CONDA_PREFIX/bin/R
```

```r
> getOption("unzip")
[1] ""
> Sys.getenv("TAR")
[1] "/bin/gtar"
> options(unzip = "/net/data.isilon/ag-saez/bq_jperales/SOFTWARE/miniconda3/envs/scmap/bin/unzip")
> Sys.setenv(TAR = "/net/data.isilon/ag-saez/bq_jperales/SOFTWARE/miniconda3/envs/scmap/bin/tar")
> devtools::install_github("mahmoudibrahim/genesorteR")
```
