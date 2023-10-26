setwd("Z:/projects/jhinckle-Love/ssGSEA")

source('common.R')
source('ssGSEAProjection.Library.R')

calculate_zscore <- function(data){
  for (i in 1:nrow(data)){
    data[i,] <- scale(data[i,])
  }
  
  return(data)
}

ssGSEA.project.dataset(input.ds = "200909Yil.gct", output.ds = "customSets_ssGSEA", gene.sets.dbfile.list = "customSets.gmx")

rc.toPlot <- calculate_zscore(toPlot)
round(mean(rowMeans(rc.toPlot)),2)
round(mean(rowVars(rc.toPlot)),2)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")