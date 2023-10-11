#Loading necessary libraries
library(WGCNA)
library(openxlsx)

options(stringsAsFactors = F)

setwd('.')


#Generate a vector of powers
powers = c(c(1:10), seq(from= 12, to= 34, by= 2))

#Power simulations
sft = pickSoftThreshold(datExpr,
                        powerVector= powers,
                        verbose= 5,
                        networkType= 'signed')

#Plotting
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

abline(h= 0.8, col= 'red')

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

sft$fitIndices

write.xlsx(sft$fitIndices, file= 'signed.xlsx')