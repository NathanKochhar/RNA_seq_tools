library(Seurat)
library(pheatmap)
library(circlize)

obj <- readRDS("/Users/nathankochhar/data/projects/1209_Phillips_Collisson/R_Analysis/main_obj.rds")


# arguments: seurat obj, group.by, ident.1, ident.2, only.pos, logfc.threshold, min.pct, recorrect_umi,
# p_val_adj.threshold, cell_type, return_matrix, cluster_rows = T, cluster_cols = T, fontsize_row = 0




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

test <- pseudobulk_DEHeatmap(obj, group.by = "condition", ident.1 = "experimental", ident.2 = "control",
                            cell.type = "cell_type")

color_mapping <- colorRamp2(c(min(test)+2, 0, max(test)-2), c("red", "white", "blue"))

pheatmap(test)

test2 <- pseudobulk_DEHeatmap(obj, group.by = "condition", ident.1 = "experimental", ident.2 = "control",
                            cell.type = "cell_type", all.genes = F, num.gene = 20)

pheatmap(test2)

test3 <- pseudobulk_DEHeatmap(obj, group.by = "condition", ident.1 = "experimental", ident.2 = "control",
                             cell.type = "cell_type", only.pos = T, logfc.threshold = 0.1, 
                             min.pct = 0.1, all.genes = F, num.gene = 20)

pheatmap(test3)

