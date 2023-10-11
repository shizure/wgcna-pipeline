library(WGCNA)
library(openxlsx)

options(stringsAsFactors = F)

setwd('.')

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

pdf(file = 'outputs/imgs/module_trait_matrix.pdf', width = 13.02, height = 7.5)

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