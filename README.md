# RNA_seq_tools
Some small tools I made because I couldn't find them online


## DEGs_per_cluster_by_condition
Small function to iteratively run FindMarkers in seurat by a chosen condition for all clusters and adds it to one data frame.

## pseudobulk_DEHeatmap
Function for a seurat object to return a pseudobulk matrix by cell type/condition for plotting with pheatmap. Celltypes will be treated as samples and log2FC is based on inputted conditions.
