library(enrichR)

# Function build off of DEGs_per_cluster_by_condition to run enrichR per cluster.
# inputs: DEGs_per_cluster_by_condition dataframe, log2FC cutoff, list of enrichR databases to use
# Will make folders and figures in working dir

PathwayAnalysis_per_cluster <- function(DEG_df, log2FC = 1, dbs) {
  
  dir.create("Pathway_Analysis")
  
  listEnrichrSites()
  setEnrichrSite("Enrichr")
  websiteLive <- TRUE
  
  DEGs_per_cell_type <- DEGs_per_cell_type[abs(DEGs_per_cell_type$avg_log2FC) > log2FC,]
  
  for (i in unique(DEGs_per_cell_type$cell_type)){
    print(i)
    geneListUp <- DEGs_per_cell_type[DEGs_per_cell_type$cell_type == i &
                                       DEGs_per_cell_type$avg_log2FC > 0,]$gene
    geneListDown <- DEGs_per_cell_type[DEGs_per_cell_type$cell_type == i &
                                         DEGs_per_cell_type$avg_log2FC < 0,]$gene
    
    #For enriched up
    if (length(geneListUp) > 5){
      dir.create(paste("Pathway_Analysis/", i, sep = ""))
      print(paste("running up regulated pathways for:", i))
      
      enriched_up <- enrichr(geneListUp, dbs)
      
      enriched_out_temp <- enriched_up
      index <- c()
      for(m in 1:length(enriched_out_temp)){
        if(nrow(enriched_out_temp[[m]]) == 0 ) {
          index <- append(index, m)}
      }
      if (is.null(index) == F){
        enriched_up <- enriched_up[-index]
      }
      
      if (length(enriched_up) != 0){
        for (j in 1:length(enriched_up)){
          title_name <- paste(names(enriched_up[j]),"_up regulated: ",i, sep = "")
          panelA<-if (websiteLive) plotEnrich(enriched_up[[j]], showTerms = 10, numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                                              title= title_name, xlab='')+
            labs( tag = "A")+theme(axis.text=element_text(size=8),plot.title= element_text(size=8,hjust=0.5))+
            theme(legend.key.size = unit(0.4, "cm"))+theme(legend.title=element_text(size=8))
          file_name <- paste("Pathway_Analysis/", i,"/up_regulated_",dbs[j],".png", sep = "")
          ggsave(file = file.path(file_name), device='png', width=7, height=4.6)
          panelA
        }
        comp_enriched_up <- bind_rows(enriched_up, .id = "enrichment_category")
        comp_enriched_up <- comp_enriched_up[order(comp_enriched_up$Adjusted.P.value),]
        write.csv(comp_enriched_up, file = paste("Pathway_Analysis/",i, "/","up_regulated_",i,"_pathways.csv", sep = ""), row.names = FALSE)
      }
    }
    
    #for enriched down
    if (length(geneListDown) > 5){
      dir.create(paste("Pathway_Analysis/", i, sep = ""))
      print(paste("running down regulated pathways for:", i))
      enriched_down <- enrichr(geneListDown, dbs)
      
      enriched_out_temp <- enriched_down
      index <- c()
      for(m in 1:length(enriched_out_temp)){
        if(nrow(enriched_out_temp[[m]]) == 0){
          index <- append(index, m)}
      }
      if (is.null(index) == F){
        enriched_down <- enriched_down[-index]
      }
      
      if (length(enriched_down) !=0){
        for (j in 1:length(enriched_down)){
          title_name <- paste(names(enriched_down[j]),"_down regulated: ",i, sep = "")
          panelA<-if (websiteLive) plotEnrich(enriched_down[[j]], showTerms = 10, numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                                              title= title_name, xlab='')+
            labs( tag = "A")+theme(axis.text=element_text(size=8),plot.title= element_text(size=8,hjust=0.5))+
            theme(legend.key.size = unit(0.4, "cm"))+theme(legend.title=element_text(size=8))
          file_name <- paste("Pathway_Analysis/", i,"/down_regulated_",dbs[j],".png", sep = "")
          ggsave(file = file.path(file_name), device='png', width=7, height=4.6)
          panelA
        }
        
        comp_enriched_down <- bind_rows(enriched_down, .id = "enrichment_category")
        comp_enriched_down <- comp_enriched_down[order(comp_enriched_down$Adjusted.P.value),]
        write.csv(comp_enriched_down, file = paste("Pathway_Analysis/",i, "/","down_regulated_",i,"_pathways.csv", sep = ""), row.names = FALSE)
      }
    }
  }
  
}

#run DEGs_per_cluster
try <- DEGs_per_cluster(obj, "condition", "experimental", "control")

#set up dbs to use
dbs_use <- c("GO_Biological_Process_2021", "KEGG_2019_Mouse","MSigDB_Hallmark_2020", "WikiPathways_2019_Mouse") 

#run function, makes folders and figures
PathwayAnalysis_per_cluster(try, dbs_use)

