## @knitr Prepare_Data

#Run a Full diagnostic for each sample
#Load the list of Seurat unfiltered objects

#Create the data list
LIST_RDS_POST_01QC_NOT_FILTERED = readRDS(PATH_01QC_OUTPUT_NOT_FILTERED)
LIST_NAMES_SAMPLES = names(LIST_RDS_POST_01QC_NOT_FILTERED)


#Loading Azimuth reference
FILE_Ref_Azimuth = file.path(PATH_EXPERIMENT_REFERENCE, "Ref_Azimuth")

#Start from the RDS of the BM object
REF_BM_OBJECT = readRDS(file.path(FILE_Ref_Azimuth, "ref_BM.Rds"))
#DimPlot(REF_BM_OBJECT, group.by = "celltype.l2")


#Build the Dataframe that will be used for NK extraction
# Create a data frame with one row and the list elements as column names
df <- as.data.frame(matrix(ncol = length(LIST_NAMES_SAMPLES), nrow = 1))
colnames(df) <- LIST_NAMES_SAMPLES

# Save the data frame as a CSV file
write.xlsx(df, file.path(PATH_EXPERIMENT_OUTPUT,"02_Azimuth_Mapping/Table_Cluster_Extraction.xlsx"), rowNames = TRUE)
