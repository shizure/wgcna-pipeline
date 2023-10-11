#Loading packages
library(WGCNA)
library(openxlsx)

source('categorical_data.R')

options(stringsAsFactors = F)

setwd('.')

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

pdf(file = 'outputs/imgs/sample_dendrogram.pdf', width = 16.66, height = 9.37)

par(cex = 0.5)
par(mar = c(6, 5, 2, 0))
plot(sampleTree, main = 'Sample clustering to detect outliers', sub = 'Average', xlab = '',
     cex.lab = 1.6, cex.axis = 1.6, cex.main = 2)
dev.off()


#Sample Dendrogram and Heatmap
traitcolors = numbers2colors(datTraits, signed = T) #Converts trait data into colors for representing in heatmap.

pdf(file = 'outputs/imgs/dendrogram_heatmap.pdf', width = 16.66, height = 9.37)

plotDendroAndColors(sampleTree, traitcolors, groupLabels = names(datTraits),
                    main = 'Sample Cluster and Heatmap', cex.lab = 1, cex.axis = 1, cex.main = 2,
                    cex.colorLabels = 0.7, cex.dendroLabels = 0.6, marAll = c(0, 6, 7, 0))

dev.off()


#Save expression data and trait data in a RData file.
save(datExpr, datTraits, file = paste(path, 'RData', 'ExprData_and_TraitsData.RData', sep = '/'))