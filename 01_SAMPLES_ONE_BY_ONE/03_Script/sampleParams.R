###############################################################################
# This file defines SAMPLE parameters as global variables that will be loaded
# before analysis starts. 
#

# Here is very basic example of the way you can declare in constants the information on samples:

# Declare all the samples names for the H5ad files
#SAMPLE_D0_NotNK = "Day0_OtherCells.h5ad"
#SAMPLE_D0_NK = "DAY0_CD56pos_NK.h5ad"
#SAMPLE_D21_IPSC_NK_DERIVED = "Day21_IPSC_NK_AFTER_EXPANSION.h5ad"

# Group the samples in a vector (usefull for loops over the samples)
#SAMPLE_H5ad_SET = c( SAMPLE_D0_NotNK, SAMPLE_D0_NK, SAMPLE_D21_IPSC_NK_DERIVED)

#SAMPLE_H5ad_SET_NoMetadata = c("Day0_OtherCells_NoMetadata.h5ad" , "DAY0_CD56pos_NK_NoMetadata.h5ad" , "Day21_IPSC_NK_AFTER_EXPANSION_NoMetadata.h5ad" )
# Define the colors that will be associated to samples all along the analysis
#SAMPLE_COLOR = c( "#00BA38", "#619CFF", "#F8766D")
#names( SAMPLE_COLOR) = SAMPLE_SET
