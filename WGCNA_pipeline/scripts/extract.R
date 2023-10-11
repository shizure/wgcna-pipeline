library(openxlsx)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(WGCNA)

options(stringsAsFactors = F)

setwd('.')

nSample = nrow(datExpr)

MEs = moduleEigengenes(datExpr, moduleColors)$eigengene

#Create trait dataframes (Can be increased if needed)
gender = as.data.frame(datTraits$Gender)
plaqueVulnerability = as.data.frame(datTraits$Plaque_Vulnerability_Index)
calsification = as.data.frame(datTraits$Calc.bin)

#Select the trait to create table
selected_trait = calsification
names(selected_trait) = strsplit(names(selected_trait), split = '\\$')[[1]][2]

#Extract module names
modNames = substring(names(MEs), 3)

#Calculates gene-module membership value and its p-value.
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = 'p'))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSample))

names(geneModuleMembership) = paste('MM', modNames, sep = '')
names(MMPvalue) = paste('p.MM', modNames, sep = '')

#Calculates gene-trait significance value and its p-value.
geneTraitSignificance = as.data.frame(cor(datExpr, selected_trait, use = 'p'))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSample))

names(geneTraitSignificance) = 'GeneSignificance'
names(GSPvalue) = 'p.GS'

#Save protein names without prefix
proteinNames = names(datExpr)
proteinNames = substring(proteinNames, 5) #5 is for ecm_ prefix, please change if needed

#Requesting gene IDs using gene names
geneIDs = mapIds(org.Hs.eg.db, keys = proteinNames,
                 column ='ENTREZID', keytype = 'SYMBOL')

#Protein names with prefix
probes = names(datExpr)

#Create temporary dataframe containing gene names, geneIDs, MM, and GS values.
geneInfo0 = data.frame(columnName = probes,
                       geneSymbol = proteinNames,
                       geneID = geneIDs,
                       module = moduleColors,
                       GS = geneTraitSignificance,
                       GSPvalue)

#Ordering modules from high absolute correlation with desired trait.
modOrder = order(-abs(cor(MEs, selected_trait, use = 'p')))

for (mod in 1:ncol(geneModuleMembership)) {
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  
  names(geneInfo0) = c(oldNames, paste('MM', modNames[modOrder[mod]], sep = ''),
                       paste('p.MM', modNames[modOrder[mod]], sep = ''))
}

#Orders the genes from high correlation to less
geneOrder = order(geneInfo0$module, -abs(geneInfo0[[names(geneTraitSignificance)]]))

#Creating final dataframe of total gene list.
geneInfo = geneInfo0[geneOrder, ]

#Save the gene info table as a xlsx file.
write.xlsx(geneInfo, 
           file = paste(output_path, 'tables', 
                        paste(names(selected_trait), '_related_gene_significance.xlsx', sep = ''), 
                        sep = '/'))