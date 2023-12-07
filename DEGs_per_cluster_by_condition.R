# Small function to iteratively run FindMarkers in seurat by a chosen condition
# for all clusters and adds it to one csv file

# The inputs are the seurat object name, the condition we're running the DE test
# over, the 2 conditions (only pairwise), whether you want only positive or both 
# positive and negitive foldchange genes, the logfc cutoff and the min.pct

DEGs_per_cluster <- function(object, condition, ident.1, ident.2, only.pos = T, 
                             logfc.threshold = 0.1, min.pct = 0.1) {
  DEGs <- data.frame()
  for (i in 0:(max(as.integer(obj@meta.data$seurat_clusters))-1)){
    subset <- subset(obj, idents = i)
    subset_degs <- FindMarkers(object = subset, group.by = condition, ident.1 = ident.1,ident.2 = ident.2,
                               only.pos = only.pos, logfc.threshold = logfc.threshold, min.pct = min.pct)
    subset_degs$gene <- rownames(subset_degs)
    subset_degs2 <- subset_degs[!grepl("^mt-", subset_degs$gene), ]
    subset_degs2 <- subset_degs[!grepl("^MT-", subset_degs$gene), ]
    subset_degs <- subset_degs[order(-abs(subset_degs$avg_log2FC)),]
    subset_degs <- subset_degs[subset_degs$p_val_adj < 0.05,]
    subset_degs$cluster <- i
    DEGs <- rbind(DEGs, subset_degs)
  }
  return(DEGs)
} 

# an example run where obj is my seurat object, "condition" is a meta.data
# variable in the object and "experimental" and "control" are the two conditions
try <- DEGs_per_cluster(obj, "condition", "experimental", "control")
try

