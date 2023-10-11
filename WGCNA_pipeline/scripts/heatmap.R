library(WGCNA)

options(stringsAsFactors = F)

setwd('.')

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