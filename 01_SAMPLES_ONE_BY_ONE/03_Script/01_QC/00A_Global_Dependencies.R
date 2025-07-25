# ##################################################
# Global declarations and libraries for the analysis
# ##################################################


######## R Libraries
library(pheatmap)
library(scales)
library(ade4)
library(Hmisc)
library(cluster)
library(ggplot2)
library(gplots)
library(gridExtra)
library(ggrepel)
library(Seurat)
library(magrittr)
library(dplyr)
library(Matrix)
library(matrixStats)
library(grid)
library(openxlsx)
library(reticulate)
library(feather)
library(reshape2)
library(plotly)
library(cowplot)
library(funr)
library(limma)
library(patchwork)
library(org.Mm.eg.db)
library(clusterProfiler)
#library(FSA)
library(Hmisc)
library(biomaRt)
library(org.Hs.eg.db)
library(biomaRt)
library(factoextra)
library(FactoMineR)
library(ggdist)
library(gghalves)
library(stringr)
library(rmdformats)
library(SingleR)
library(celldex)
library(corrplot)
library(RColorBrewer)
library(readxl)
library(rlist)



#Def functions manualy
#def function
# +- std
data_summary <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}





