library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)

options(stringsAsFactors = F)

setwd('.')


#Modules that are used in enrichment.
modules = substring(names(MEs), 3) #If all modules are not needed, please create a list of desired module names.

# Example of module list
# modules = c('pink', 'blue', 'black', 'yellow',
#             'red', 'green', 'magenta', 'turquoise', 'brown')


#Create module folders to save outputs
for(module in modules){
  if(!dir.exists(paste('outputs/tables', module, sep = '/'))){
    dir.create(paste('outputs/tables', module, sep = '/'))
  }
}


#GO enrichment for each module above
for(i in modules){
  go_enrichment = enrichGO(gene = geneInfo[geneInfo$module == i, ]$geneID, 
                           OrgDb ='org.Hs.eg.db',
                           keyType = "ENTREZID",
                           ont = "BP",  # Biological Process ontology
                           pAdjustMethod = "BH",  # Adjust p-values
                           pvalueCutoff = 0.05,
                           readable = T)
  
  result = data.frame(GO_ID = go_enrichment$ID,
                      GO_Description = go_enrichment$Description,
                      P.value = go_enrichment$pvalue,
                      GeneID = go_enrichment$geneID)
  
  write.xlsx(result[order(result$P.value), ],
             file = paste('outputs/tables/', i, '/GOenrichment.xlsx', sep = '/'))
  
  
  print(paste('Generating ', i, ' module GO enrichment table is completed.'))
  
  
  
  
  #KEGG pathway enrichment for each module above
  kegg_enrichment = enrichKEGG(gene = geneInfo[geneInfo$module == i, ]$geneID,
                               organism = 'hsa',
                               keyType = 'ncbi-geneid',
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05)
  
  result = data.frame(KEGG_ID = kegg_enrichment$ID,
                      KEGG_Description = kegg_enrichment$Description,
                      P.value = kegg_enrichment$pvalue,
                      GeneID = kegg_enrichment$geneID)
  
  write.xlsx(result[order(result$P.value), ],
             file = paste('outputs/tables/', i, '/KEGGenrichment.xlsx', sep = '/'))
  
  print(paste('Generating ', i, ' module KEGG enrichment table is completed.'))
  
  
  
  #Gene List
  write.xlsx(geneInfo[geneInfo$module == i, c(1:4)],
             file = paste('outputs/tables/', i, '/gene_list.xslx', sep = '/'))
}