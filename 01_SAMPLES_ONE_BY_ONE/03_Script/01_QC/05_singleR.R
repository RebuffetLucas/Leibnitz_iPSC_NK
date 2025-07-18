## @knitr SingleR

# ---- Section header for RMarkdown output ----
cat(" \n \n")
cat("#### Automatic identification with SingleR ")
cat(" \n \n")

# ---- Convert Seurat object to SingleCellExperiment format ----
# Required format for SingleR input
sce <- as.SingleCellExperiment(myObjectSeurat)

# ---- Load reference data ----
# BlueprintEncodeData contains immune cell types (fine-level resolution)
BlueprintEncodeData.ref <- celldex::BlueprintEncodeData()

# ---- Run SingleR to annotate cells ----
# Matches expression profiles from 'sce' to the reference set using fine-level labels
BlueprintEncodeData.fine <- SingleR(
  test = sce,
  ref = BlueprintEncodeData.ref,
  labels = BlueprintEncodeData.ref$label.fine
)

# ---- Store SingleR results in Seurat metadata ----
myObjectSeurat@meta.data$BlueprintEncodeData.fine <- BlueprintEncodeData.fine$labels

# ---- Prepare summary table of cell annotations ----
Table_Automatic_Id <- data.frame(table(myObjectSeurat$BlueprintEncodeData.fine))
colnames(Table_Automatic_Id) <- c("SingleR_Id", "Number_of_Cells")

# ---- Define color palette ----
Nb_markers <- length(unique(Table_Automatic_Id$SingleR_Id))
mypalette <- hue_pal()(Nb_markers)

# ---- Name colors by cell type ----
names(mypalette) <- levels(Table_Automatic_Id$SingleR_Id)

# ---- Highlight NK cells in red ----
mypalette["NK cells"] <- "#FF0000"

# ---- Plot UMAP colored by SingleR-assigned identity ----
cat("#####" ,"UMAP","\n")

cat(' \n \n')
print(
  DimPlot(
    myObjectSeurat,
    group.by = "BlueprintEncodeData.fine",
    label = TRUE,
    label.size = LABEL_SIZE,
    cols = mypalette
  ) &
    ggtitle("Automatic identification with SingleR")
)
cat(' \n \n')


# ---- Display annotation table as interactive DataTable ----
cat(' \n \n')
print(
  htmltools::tagList(
    DT::datatable(
      Table_Automatic_Id,
      rownames = FALSE,
      extensions = 'Buttons',
      options = list(
        dom = 'Blfrtip',
        buttons = c('excel', "csv"),
        fixedHeader = TRUE
      )
    ) %>%
      DT::formatStyle(
        'SingleR_Id',
        backgroundColor = DT::styleEqual(
          sort(unique(Table_Automatic_Id$SingleR_Id)),
          mypalette[1:Nb_markers]
        )
      )
  )
)
cat(' \n \n')

# ---- Barplot: SingleR cell types per Seurat cluster ----
# Stacked barplot showing composition of each Seurat cluster by cell type
p5 <- ggplot(myObjectSeurat@meta.data, aes(x = seurat_clusters, fill = BlueprintEncodeData.fine)) +
  geom_bar(position = "fill") +
  theme(
    text = element_text(size = 20),
    axis.text.x = element_text(angle = 90)
  ) +
  scale_fill_manual(values = mypalette[1:Nb_markers])

cat(' \n \n')
print(p5)
cat(' \n \n')




## @knitr CellsToRemove

# ---- Identify cells to remove based on SingleR labels ----
cat("#####" ,"Highlight cells to remove","\n")

# Set cell identities to SingleR fine annotations
myObjectSeurat <- SetIdent(myObjectSeurat, value = "BlueprintEncodeData.fine")

# Select cells from unwanted lineages (e.g. myeloid/monocytic contamination)
cells <- WhichCells(myObjectSeurat, idents = c("MEP", "GMP", "Monocytes"))

# Highlight these cells on the UMAP
print(
  DimPlot(
    myObjectSeurat,
    reduction = "umap",
    cells.highlight = cells,
    pt.size = UMAP_PT_SIZE,
    sizes.highlight = 1.5
  ) +
    scale_color_manual(
      labels = c("others", "cells to remove"),
      values = c("grey", "#CE4C4B")  # Red for cells to remove
    )
)
cat(' \n \n')


# ---- Display barcodes of cells to remove ----
cat("#####" ,"Barcodes","\n")

# Print interactive table of cell barcodes for flagged cells
DT::datatable(
  as.data.frame(cells),
  rownames = FALSE,
  extensions = 'Buttons',
  options = list(
    dom = 'Blfrtip',
    buttons = c('excel', "csv"),
    fixedHeader = TRUE
  )
)
cat(' \n \n')

