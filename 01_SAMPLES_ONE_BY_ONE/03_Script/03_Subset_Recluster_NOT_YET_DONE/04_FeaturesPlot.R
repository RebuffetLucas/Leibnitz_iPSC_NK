## @knitr FeaturePlot

cat(" \n \n")
cat("#### FeaturePlot ")
cat(" \n \n")

cat(" \n \n")
print(  FeaturePlot(myObjectSeurat,features = c("percent.mito","percent.ribo","nFeature_RNA","nCount_RNA"),pt.size=1.2)& scale_color_viridis_c(option="plasma") )
cat(" \n \n")


cat(" \n \n")
print(  FeaturePlot(myObjectSeurat,features = c("FCER1G" ,"SPON2" , "GZMB",
                                                "XCL1","XCL2", "GZMK",
                                                "CD3E","CD3D","KLRC2" ),pt.size=1.2)& scale_color_viridis_c(option="plasma") )
cat(" \n \n")


cat(" \n \n")
print(  FeaturePlot(myObjectSeurat,features = c("NKG7","LYZ", "GZMB",
                                                "XCL1","XCL2", "GZMK",
                                                "CD3E","CD3D","KLRC2" ),pt.size=1.2)& scale_color_viridis_c(option="plasma") )
cat(" \n \n")






## @knitr Highlight_cells_to_remove

cat("#####" ,"UMAP","\n")

cells <- WhichCells(myObjectSeurat, idents=3)
print(DimPlot(myObjectSeurat,reduction = "umap",cells.highlight=cells,pt.size = 1,sizes.highlight = 1.5) + scale_color_manual(labels = c("others", "cells to remove"), values = c("grey","#CE4C4B")))
cat(' \n \n')

cat("#####" ,"Barcodes","\n")

DT::datatable(as.data.frame(cells) ,rownames = FALSE,extensions = 'Buttons', 
              options = list(dom = 'Blfrtip', 
                             buttons = c('excel', "csv"), fixedHeader = TRUE)
)
cat(' \n \n')





## @knitr Scoring_NK123_and_13Genes


cat(" \n \n")
cat("#### Scoring NK 1, NK2, NK3 and 13 Genes Crinier ")
cat(" \n \n")


  #NK1 , NK2 and NK3
Markers_Seurat= readRDS( file.path(PATH_EXPERIMENT_REFERENCE ,"Tables_Gene_Signatures/All_Markers_MetaNK_CITEseq3clusters.rds" )) #FindAllMarkersOutput format


#Extracting top markers
Markers_Seurat %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but not for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC) -> top_All


#Building list of best markers
list_Top_Genes = list()
for (i in levels(top_All$cluster)){
  top_All %>%
    filter(cluster == i ) -> top_clust
  list_Top_Genes  =append(list_Top_Genes, list(top_clust$gene))
}

names(list_Top_Genes) = levels(top_All$cluster)
#print(list_Top_Genes)


#Format that can be used by AddModuleScore
List_To_Use = lapply(list_Top_Genes , function(x) as.data.frame(x))
MONITORED_Markers = List_To_Use

for (i in names(MONITORED_Markers)){
  myObjectSeurat = AddModuleScore(myObjectSeurat, features = as.list(MONITORED_Markers[[i]]), pool= NULL ,name= i , seed=19)
}

#VlnPlot and save

p2 = DimPlot(myObjectSeurat, reduction = "umap")

p3 = FeaturePlot(myObjectSeurat, features = paste0(names(MONITORED_Markers),"1"), max.cutoff = 1.5)   &    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
p4 = VlnPlot(myObjectSeurat, features = paste0(names(MONITORED_Markers),"1"), y.max = 1.5, pt.size = 0, ncol = 1) & stat_summary(fun.data=data_summary,color="black")

cat(' \n \n')
print(p3 + p2)
cat(' \n \n')

cat(' \n \n')
print(p4)
cat(' \n \n')



#Crinier 2018 13 genes

FILE_SIGNATURES_PATH = file.path(PATH_EXPERIMENT_REFERENCE ,"Tables_Gene_Signatures","Crinier_13_NKgenes.xlsx" ) #FindAllMarkersOutput format


###Markers for scoring of subsets
# getting data from sheets
sheets_names <- openxlsx::getSheetNames(FILE_SIGNATURES_PATH)

#def function

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}


MONITORED_Markers = read_excel_allsheets(FILE_SIGNATURES_PATH)
MONITORED_Markers = lapply(MONITORED_Markers, FUN=function(x) intersect(unlist(x), rownames(myObjectSeurat)))
MONITORED_Markers= lapply(MONITORED_Markers, FUN=function(x) head(x, n=20))

myObjectSeurat = AddModuleScore(myObjectSeurat, features = MONITORED_Markers, name= "Crinier_13genes_Score" , pool= NULL , seed=19)
p5 = VlnPlot(myObjectSeurat, features = "Crinier_13genes_Score1",  y.max = 1.5, pt.size = 0 ) & stat_summary(fun.data=data_summary,color="black")

cat(' \n \n')
print(p5)
cat(' \n \n')

