## @knitr SingleR

cat(" \n \n")
cat("#### Automatic identification with SingleR ")
cat(" \n \n")


sce <- as.SingleCellExperiment(myObjectSeurat)
BlueprintEncodeData.ref <- celldex::BlueprintEncodeData()
BlueprintEncodeData.fine <- SingleR(test = sce,ref = BlueprintEncodeData.ref,labels = BlueprintEncodeData.ref$label.fine)

myObjectSeurat@meta.data$BlueprintEncodeData.fine <- BlueprintEncodeData.fine$labels



#Prepare table
#print tables of automatic identification with singleR

Table_Automatic_Id <- data.frame( table(myObjectSeurat$BlueprintEncodeData.fine))
colnames(Table_Automatic_Id) = c("SingleR_Id" , "Number_of_Cells")
Nb_markers=length(unique(Table_Automatic_Id$SingleR_Id))
mypalette= hue_pal()(Nb_markers)

#Turn NK cells to red
# Create a named vector for the palette
names(mypalette) <- levels(Table_Automatic_Id$SingleR_Id)
# Change the color of "NK cells" to #FF0000
mypalette["NK cells"] <- "#FF0000"


#Print the UMAP
cat("#####" ,"UMAP","\n")

cat(' \n \n')
print(DimPlot(myObjectSeurat,group.by = "BlueprintEncodeData.fine", label = TRUE, label.size = LABEL_SIZE, cols = mypalette) & ggtitle("Automatic identification with SingleR"))
cat(' \n \n')




#Print the table
cat(' \n \n')
print( htmltools::tagList(DT::datatable(Table_Automatic_Id,rownames = FALSE,extensions = 'Buttons', 
                                        options = list(dom = 'Blfrtip', 
                                                       buttons = c('excel', "csv"), fixedHeader = TRUE)
)%>% 
  DT::formatStyle(
    'SingleR_Id',
    backgroundColor = DT::styleEqual(sort(unique(Table_Automatic_Id$SingleR_Id)),  mypalette[1:Nb_markers])
  )))

cat(' \n \n')


#print barplot of cluster composition

p5 = ggplot(myObjectSeurat@meta.data, aes(x=seurat_clusters, fill= BlueprintEncodeData.fine)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= mypalette[1:Nb_markers])

cat(' \n \n')
print(p5)
cat(' \n \n')






## @knitr CellsToRemove
cat("#####" ,"Highlight cells to remove","\n")

myObjectSeurat <- SetIdent(myObjectSeurat, value = "BlueprintEncodeData.fine")

cells <- WhichCells(myObjectSeurat, idents=c("MEP","GMP","Monocytes"))
print(DimPlot(myObjectSeurat,reduction = "umap",cells.highlight=cells,pt.size = 1,sizes.highlight = 1.5) + scale_color_manual(labels = c("others", "cells to remove"), values = c("grey","#CE4C4B")))
cat(' \n \n')


cat("#####" ,"Barcodes","\n")

DT::datatable(as.data.frame(cells) ,rownames = FALSE,extensions = 'Buttons', 
              options = list(dom = 'Blfrtip', 
                             buttons = c('excel', "csv"), fixedHeader = TRUE)
)
cat(' \n \n')
