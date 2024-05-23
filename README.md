# scRNA-seq Tools
Some small tools I made because I couldn't find them online, mostly seurat based

#

### DEGs_per_cluster_by_condition
Small function to iteratively run FindMarkers in seurat by a chosen condition for all clusters and adds it to one data frame. Can easily be modified to find DEGs per specified meta.data column (cell type).
```R
library(Seurat)
# an example run where obj is my seurat object, "condition" is a meta.data
# variable in the object and "experimental" and "control" are the two conditions

try <- DEGs_per_cluster(obj, "condition", "experimental", "control")
try
```
#
### pseudobulk_DEHeatmap
Function for a seurat object to return a pseudobulk matrix by cell type/condition for plotting with pheatmap. Celltypes will be treated as samples and log2FC is based on inputted conditions.
```R
library(Seurat)
library(pheatmap)
library(circlize)
# an example run where obj is my seurat object, "condition" is a meta.data
# variable in the object and "experimental" and "control" are the two conditions
# and cell type is "cell_type"

test <- pseudobulk_DEHeatmap(obj, group.by = "condition", ident.1 = "experimental", ident.2 = "control",
                            cell.type = "cell_type")

pheatmap(test)
```
#
### PathwayAnalysis_per_cluster
Function to preform EnrichR analysis per cluster based off of the DEGs_per_cluster_by_condition function output. Will create a folder for each cluster and will make EnrichR barplots per dbs per cluster. Can be editted along with DEGs_per_cluster_by_condition to run pathway analysis based on specified meta.data column (cell type).
```R
library(enrichR)

#run DEGs_per_cluster
try <- DEGs_per_cluster(obj, "condition", "experimental", "control")

#set up dbs to use
dbs_use <- c("GO_Biological_Process_2021", "KEGG_2019_Mouse","MSigDB_Hallmark_2020", "WikiPathways_2019_Mouse") 

#run function, makes folders and figures
PathwayAnalysis_per_cluster(try, dbs_use)
```
