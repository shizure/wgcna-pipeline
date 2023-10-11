#Loading necessary libraries
library(WGCNA)
library(openxlsx)

options(stringsAsFactors = F)

setwd('.')


#Network Construction Using blockwiseModules() function
power = 14
minModuleSize = 10
seed = 12345

net = blockwiseModules(datExpr,
                       maxBlockSize = 10000,
                       randomSeed = seed,
                       minModuleSize = minModuleSize,
                       power = power,
                       networkType = 'signed',
                       TOMType = 'signed',
                       deepSplit = 3,
                       mergeCutHeight = 0.3,
                       nThreads = 10,
                       verbose = 3,
                       saveTOMs = F)

table(net$colors)

#Plotting Modules
mergedColors = labels2colors(net$colors)

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    'Module Colors', dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05,
                    main = 'Module Cluster')

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
proteinTree = net$dendrograms[[1]]


save(datExpr, datTraits, net, seed, minModuleSize, moduleColors, moduleLabels, proteinTree,
     file= 'pwr14.RData')