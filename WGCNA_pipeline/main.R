#####NETWORK CONSTRUCTION AND MODULE ANALYSIS PIPELINE#####
#To run the script from terminal, add Rscript to your PATH.
#Pipeline requires only two arguments; working directory and directory of the dataset
#Dataset must contain clinical trait and expression data together
#Change row 61 according to your dataset
#See categorical_data.R to adjust categorical values in your dataset
#Most part of the pipeline have own plot functions to avoid repeat. 
#If any change in plot function is needed, see the plot functions in associated part.
#The pipeline can be used several times in a row, the saved files will be named with unique number suffixes.
#To avoid any directory confussion, please keep your data in 'data' folder, and please keep your data folder in the same level with the scripts.
#!!!!!Before re-run the pipeline, please save outputs/tables directory. It will be overwritten.
########### WORKSPACE SETTINGS ##########
setwd('.')
#Load packages
library(openxlsx)
library(org.Hs.eg.db)
library(clusterProfiler)
library(WGCNA)

#Load external script
source('categorical_data.R')

#Set workspace pathway
path = getwd()

options(stringsAsFactors = F)

#Check if necessary folders exist in path; if not create
if(!dir.exists(paste(path, 'outputs', sep = '/'))){
  dir.create(paste(path, 'outputs', sep = '/'))
}
if(!dir.exists(paste(path, 'outputs/imgs', sep = '/'))){
  dir.create(paste(path, 'outputs/imgs', sep = '/'))
}
if(!dir.exists(paste(path, 'outputs/tables', sep = '/'))){
  dir.create(paste(path, 'outputs/tables', sep = '/'))
}
if(!dir.exists(paste(path, 'RData', sep = '/'))){
  dir.create(paste(path, 'RData', sep = '/'))
}
output_path = paste(path, 'outputs', sep = '/')


#Create a resolution matrix. Plotting functions use this matrix.
resolutions = data.frame(width = c(13.02, 16.66),
                         height = c(7.5, 9.37))
rownames(resolutions) = c('small', 'large')





########## *1* PREPROCESSING #########

#Read the dataset containing expression and trait data
dataset = read.xlsx('data/MergedPatientData.xlsx')


##Seperate Expression data and Trait data
#Expression data#
datExpr = as.data.frame(dataset[, -c(1:49)]) #use t() function to make sure proteins are in columns.
rownames(datExpr) = dataset$STUDY_NUMBER     #Attach patiend IDs to the rows.


#Trait data#
allTraits = categorical_to_numeric(dataset) #Please see categorical_data.R script to modify.

#Check each column if it is empty or not; if so removes the column and prints which columns were removed.
for(i in names(allTraits)){
  if (all(is.na(allTraits[[i]]))){
    allTraits[[i]] = NULL
    print(sprintf('%s was an empty feature. Thus, it was removed from traits dataset.', i))
  }
}

#Match and sort expression and trait datasets
patientSamples = rownames(datExpr)
traitRows = match(patientSamples, allTraits$STUDY_NUMBER)
datTraits = allTraits[traitRows, -1]
collectGarbage()

#Use for trait elimination.
#datTraits = datTraits[, -c(3, 9, 14, 15, (17:20), (21:23), 27, 29, 30, 42)]

#Check for Null values
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK #if returns true continue; if not, you should handle missing values.


##Generate sample dendrogram and heatmap figures.
#Hierarchical clustering of samples to detect outliers.
sampleTree = hclust(dist(datExpr), method = 'average')


#PLOT FUNCTION#
plot_sampleDendrogram = function(filename, res) {
  pdf(file = paste(output_path, 'imgs', filename, sep = '/'), 
      width = resolutions[res, 1], height = resolutions[res, 2], )
  
  par(cex = 0.5)
  par(mar = c(6, 5, 2, 0))
  plot(sampleTree, main = 'Sample clustering to detect outliers', sub = 'Average', xlab = '',
       cex.lab = 1.6, cex.axis = 1.6, cex.main = 2)
  dev.off()
}

if(!file.exists(paste(output_path, 'imgs', 'sample_dendrogram.pdf', sep = '/'))){
  plot_sampleDendrogram('sample_dendrogram.pdf', 'large')
} else {
  suffix = 1
  repeat{
    new_name = paste('sample_dendrogram(', suffix, ').pdf', sep = '')
    if(!file.exists(paste(output_path, 'imgs', new_name, sep = '/'))){
      plot_sampleDendrogram(new_name, 'large')
      break
    } else{
      suffix = suffix + 1
    }
  }
}


#Sample Dendrogram and Heatmap
traitcolors = numbers2colors(datTraits, signed = T) #Converts trait data into colors for representing in heatmap.

#PLOT FUNCTION#
plot_dendrogramHeatmap = function(filename, res){
  pdf(file = paste(output_path, 'imgs', filename, sep = '/'), 
       width = resolutions[res, 1], height = resolutions[res, 2])
  plotDendroAndColors(sampleTree, traitcolors, groupLabels = names(datTraits),
                      main = 'Sample Cluster and Heatmap', cex.lab = 1, cex.axis = 1, cex.main = 2,
                      cex.colorLabels = 0.7, cex.dendroLabels = 0.6, marAll = c(0, 6, 7, 0))
  dev.off()
}

if(!file.exists(paste(output_path, 'imgs', 'dendrogram_heatmap.pdf', sep = '/'))){
  plot_dendrogramHeatmap('dendrogram_heatmap.pdf', 'large')
} else {
  suffix = 1
  repeat{
    new_name = paste('dendrogram_heatmap(', suffix, ').pdf', sep = '')
    if(!file.exists(paste(output_path, 'imgs', new_name, sep = '/'))){
      plot_dendrogramHeatmap(new_name, 'large')
      break
    } else{
      suffix = suffix + 1
    }
  }
}



##Save expression data and trait data in a RData file.
#List of objects that goint to be saved
save_list = c(datExpr, datTraits)

#Saves the save_list in an unique named RData file
if(!file.exists(paste(path, 'RData', 'ExprData_and_TraitsData.RData', sep = '/'))){
  save(save_list, file = paste(path, 'RData', 'ExprData_and_TraitsData.RData', sep = '/'))
} else {
  suffix = 1
  repeat{
    new_name = paste('ExprData_and_TraitsData(', suffix, ').RData', sep = '')
    if(!file.exists(paste(path, 'RData', new_name, sep = '/'))){
      save(save_list, file = paste(path, 'RData', new_name, sep = '/'))
      break
    } else{
      suffix = suffix + 1
    }
  }
}






########## *2* SOFT-THRESHOLD (POWER) ANALYSIS ##########

#Create a vector of powers
powers = c(c(1:10), seq(from= 12, to= 34, by= 2))

#Power simulations
sft = pickSoftThreshold(datExpr,
                        powerVector= powers,
                        verbose= 5,
                        networkType= 'signed')



#Plotting
#PLOT FUNCTION#
plot_softThreshold = function(filename, res){
  pdf(file = paste(output_path, 'imgs', filename, sep = '/'), 
       width = resolutions[res, 1], height = resolutions[res, 2])
  
  par(mfrow= c(1, 2))
  
  plot(sft$fitIndices[, 1],
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab= 'SoftThreshold(power)',
       ylab= 'ScaleFreeTopologyModelFit,signedR^2',
       type= 'n', main= 'Scaleindependence')
  
  text(sft$fitIndices[,1],
       -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
       labels= powers
       ,col= 'black',
       cex= 0.8)
  
  abline(h= 0.8, col= 'red') #to identify the threshold (not necessary)
  
  plot(sft$fitIndices[,1],
       sft$fitIndices[,5], 
       xlab= 'SoftThreshold(power)',
       ylab= 'MeanConnectivity',
       type= 'n', main= 'Meanconnectivity')
  
  text(sft$fitIndices[,1],
       sft$fitIndices[,5],
       labels= powers,
       col='red',
       cex= 0.8)
  dev.off()
}

if(!file.exists(paste(output_path, 'imgs', 'soft_threshold.pdf', sep = '/'))){
  plot_softThreshold('soft_threshold.pdf', 'small')
} else {
  suffix = 1
  repeat{
    new_name = paste('soft_threshold(', suffix, ').pdf', sep = '')
    if(!file.exists(paste(output_path, 'imgs', new_name, sep = '/'))){
      plot_softThreshold(new_name, 'small')
      break
    } else{
      suffix = suffix + 1
    }
  }
}


#If wanted, sft$fitIndices might be saved as a xlsx file using write.xlsx().






########## *3* NETWORK CONSTRUCTION #########

#Network Construction Using blockwiseModules() function
power = 14 #The value selected from soft-threshold analysis.
minModuleSize = 10 #Minimum number of genes in each module.
seed = 12345 #For repeatability. 

#There are more parameters in blockWiseModules(). Please see them before usage.
net = blockwiseModules(datExpr,                       #Expression data.
                       maxBlockSize = 10000,          #Maximum size of blocks that are processed at onetime.
                       randomSeed = seed,
                       minModuleSize = minModuleSize,
                       power = power,
                       networkType = 'signed',        #'signed', 'unsigned', or 'signed hybrid'
                       TOMType = 'signed',
                       deepSplit = 3,                 #Sensitivity between value 0 to 4.
                       mergeCutHeight = 0.3,          #Threshold where modules below are merged.
                       nThreads = 10,                 #Number of allowed threads from CPU. 
                       verbose = 3)                   #Verbosity level between 0 to 7


##Plotting Modules##
moduleColors = labels2colors(net$colors)


#PLOT FUNCTION#
plot_moduleDendrogram = function(filename, res){
  pdf(file = paste(output_path, 'imgs', filename, sep = '/'), 
       width = resolutions[res, 1], height = resolutions[res, 2])
  
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      'Module Colors', dendroLabels = F, hang = 0.03,
                      addGuide = T, guideHang = 0.05,
                      main = 'Module Cluster')
  dev.off()
}

if(!file.exists(paste(output_path, 'imgs', 'module_dendrogram.pdf', sep = '/'))){
  plot_moduleDendrogram('module_dendrogram.pdf', 'large')
} else {
  suffix = 1
  repeat{
    new_name = paste('module_dendrogram(', suffix, ').pdf', sep = '')
    if(!file.exists(paste(output_path, 'imgs', new_name, sep = '/'))){
      plot_moduleDendrogram(new_name, 'large')
      break
    } else{
      suffix = suffix + 1
    }
  }
}



#Variables will be used further steps.
moduleLabels = net$colors
MEs = moduleEigengenes(datExpr, moduleColors, softPower = power)$eigengenes
proteinTree = net$dendrograms[[1]]


##Save important variables in a RData file.
#List of objects that goint to be saved
save_list = c(datExpr, datTraits, net, seed, minModuleSize, moduleColors, moduleLabels, MEs, proteinTree)

#Saves the save_list in an unique named RData file
if(!file.exists(paste(path, 'RData', 'network.RData', sep = '/'))){
  save(save_list, file = paste(path, 'RData', 'network.RData', sep = '/'))
} else {
  suffix = 1
  repeat{
    new_name = paste('network.RData(', suffix, ').RData', sep = '')
    if(!file.exists(paste(path, 'RData', new_name, sep = '/'))){
      save(save_list, file = paste(path, 'RData', new_name, sep = '/'))
      break
    } else{
      suffix = suffix + 1
    }
  }
}





########## *4* MODULE-TRAIT MATRIX #########

nSample = nrow(datExpr)

#Order MEs
MEs = orderMEs(MEs) #It also can be used to see Module-Sample correlation.

#Calculates module and trait correlation using Pearson correlation.
moduleTraitCor = cor(MEs, datTraits, use = 'p')
#Calculates p-value of correlation
moduleTraitPValue = corPvalueStudent(moduleTraitCor, nSample) 


##Plotting##

textMatrix = paste(signif(moduleTraitCor, 2), '\n(',
                   signif(moduleTraitPValue, 1), ')', sep = '')

dim(textMatrix) = dim(moduleTraitCor)


#PLOT FUNCTION#
plot_matrix = function(filename, res){
  pdf(file = paste(output_path, 'imgs', filename, sep = '/'), 
       width = resolutions[res, 1], height = resolutions[res, 2])
  
  par(mar = c(6, 8.5, 3, 3))
  
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = F,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = F,
                 cex.text = 0.5,
                 cex.lab.x = 0.7,
                 zlim = c(-1, 1),
                 main = 'Module-Trait Relationship')
  dev.off()
}

if(!file.exists(paste(output_path, 'imgs', 'module_trait_matrix.pdf', sep = '/'))){
  plot_matrix('module_trait_matrix.pdf', 'small')
} else {
  suffix = 1
  repeat{
    new_name = paste('module_trait_matrix(', suffix, ').pdf', sep = '')
    if(!file.exists(paste(output_path, 'imgs', new_name, sep = '/'))){
      plot_matrix(new_name, 'small')
      break
    } else{
      suffix = suffix + 1
    }
  }
}





########## *5* MODULE HEATMAP AND EIGENGENE DENDROGRAM ##########

#Generate module heatmap using TOM dissimilarity
dissTOM = 1 - TOMsimilarityFromExpr(datExpr, power = 14, networkType = 'signed', verbose = 3)

plotTOM = dissTOM^7 #Increase the power to clear visual
diag(plotTOM) = NA #Marks diagonal to better visualization

#PLOT FUNCTION#

plot_networkHeatmap = function(filename, res){
  pdf(file = paste(output_path, 'imgs', filename, sep = '/'), 
       width = resolutions[res, 1], height = resolutions[res, 2])
  
  TOMplot(plotTOM, proteinTree, moduleColors, main = 'Network Heatmap plot, all genes')
  dev.off()
}

if(!file.exists(paste(output_path, 'imgs', 'network_heatmap.pdf', sep = '/'))){
  plot_networkHeatmap('network_heatmap.pdf', 'small')
} else {
  suffix = 1
  repeat{
    new_name = paste('network_heatmap(', suffix, ').pdf', sep = '')
    if(!file.exists(paste(output_path, 'imgs', new_name, sep = '/'))){
      plot_networkHeatmap(new_name, 'small')
      break
    } else{
      suffix = suffix + 1
    }
  }
}


#Eigengene Dendrogram and Module-Module Correlation

#PLOT FUNCTION#
plot_eigengeneNetwork = function(filename, res){
  pdf(file = paste(output_path, 'imgs', filename, sep = '/'), 
       width = resolutions[res, 1], height = resolutions[res, 2])
  
  plotEigengeneNetworks(MEs, 'Eigengene Dendrogram', 
                        plotHeatmaps = F, marDendro = c(0, 4, 1, 2), cex.lab = 1)
  dev.off()
}

if(!file.exists(paste(output_path, 'imgs', 'eigengene_dendrogram.pdf', sep = '/'))){
  plot_eigengeneNetwork('eigengene_network.pdf', 'small')
} else {
  suffix = 1
  repeat{
    new_name = paste('eigengene_dendrogram(', suffix, ').pdf', sep = '')
    if(!file.exists(paste(output_path, 'imgs', new_name, sep = '/'))){
      plot_eigengeneNetwork(new_name, 'small')
      break
    } else{
      suffix = suffix + 1
    }
  }
}

#PLOT FUNCTION#
plot_eigengeneAdjacency = function(filename, res){
  pdf(file = paste(output_path, 'imgs', filename, sep = '/'), 
       width = resolutions[res, 1], height = resolutions[res, 2])
  
  plotEigengeneNetworks(MEs, 'Eigengene Adjacency',
                        plotDendrograms = F, marHeatmap = c(3, 4, 1, 2), cex.lab = 1, 
                        signed = T, heatmapColors = greenWhiteRed(50), xLabelsAngle = 90)
  dev.off()
}

if(!file.exists(paste(output_path, 'imgs', 'eigengene_adjacency.pdf', sep = '/'))){
  plot_eigengeneAdjacency('eigengene_adjacency.pdf', 'small')
} else {
  suffix = 1
  repeat{
    new_name = paste('eigengene_adjacency(', suffix, ').pdf', sep = '')
    if(!file.exists(paste(output_path, 'imgs', new_name, sep = '/'))){
      plot_eigengeneAdjacency(new_name, 'small')
      break
    } else{
      suffix = suffix + 1
    }
  }
}





########## *6* EXTRACTING GENES ##########

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



# ##Use part of script below to generate MM vs GS plot for a desired module.
# module = 'brown' #Module that will be used to plot
# column = match(module, modNames) #Selects module genes from all genes.
# moduleGenes = moduleColors == module
# 
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste('Module membership in', module, 'module'),
#                    ylab = 'Gene significance for plaque vulnerability',
#                    main = paste('Module membership vs Gene significance\n'),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# 
# #Lists the gene names of desired module.
# names(datExpr)[moduleColors == 'brown']





########## *7* GO AND KEGG ENRICHMENT #########

#Modules that are used in enrichment.
modules = substring(names(MEs), 3) #If all modules are not needed, please create a list of desired module names.

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
             file = paste(output_path, 'tables', i, 'GOenrichment.xlsx', sep = '/'))
  
  
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
             file = paste(output_path, 'tables', i, 'KEGGenrichment.xlsx', sep = '/'))
  
  print(paste('Generating ', i, ' module KEGG enrichment table is completed.'))
  
  
  
  #Gene List
  write.xlsx(geneInfo[geneInfo$module == i, c(1:4)],
             file = paste(output_path, 'tables', i, 'gene_list.xlsx', sep = '/'))
}