# arguments: seurat obj, group.by, ident.1, ident.2, only.pos, logfc.threshold, min.pct, recorrect_umi,
# p_val_adj.threshold, cell_type, all.genes, num.gene

pseudobulk_DEHeatmap <- function(obj, group.by, ident.1, ident.2, cell.type, only.pos = FALSE, logfc.threshold = 1, 
                                min.pct = 0.05, recorrect_umi = FALSE, all.genes = T, num.gene = 0){
  table(obj$orig.ident)
  bulk_degs <- FindMarkers(object = obj, group.by = group.by, ident.1 = ident.1, ident.2 = ident.2, only.pos = only.pos, 
                           logfc.threshold = logfc.threshold, min.pct = min.pct, recorrect_umi = recorrect_umi)
  bulk_degs$gene <- rownames(bulk_degs)
  bulk_degs <- bulk_degs[bulk_degs$p_val_adj < 0.05,]
  bulk_degs <- bulk_degs[!grepl("^mt-", bulk_degs$gene), ]
  bulk_degs <- bulk_degs[!grepl("^MT-", bulk_degs$gene), ]
  
  bulk_degs <- bulk_degs[order(abs(bulk_degs$avg_log2FC), decreasing = TRUE), ]
  
  if(all.genes == T){
    genes = bulk_degs$gene
  } else {
    genes = bulk_degs$gene[1:num.gene]
  }
  
  mat <- data.frame(row.names = genes)
  
  obj@meta.data$cell_type1 = obj@meta.data[[cell.type]]
  
  for (i in levels(obj@meta.data[[cell.type]])){
    print(i)
    subset <- subset(obj, subset = cell_type1 == i)
    subset_degs <- FindMarkers(object = subset, group.by = group.by, ident.1 = ident.1, ident.2 = ident.2,
                               feature = genes, logfc.threshold = 0, min.pct = 0, recorrect_umi = FALSE)
    subset_degs = subset_degs[,2, drop = FALSE]
    colnames(subset_degs) = i
    mat <- merge(mat, subset_degs, by = 'row.names')
    rownames(mat) <- mat$Row.names
    mat <- mat[,-1, drop = FALSE]
  }
  
  mat <- data.matrix(mat)
  
  return (mat)
}
