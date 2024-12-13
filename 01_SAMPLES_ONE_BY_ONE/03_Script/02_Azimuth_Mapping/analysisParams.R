###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#

ANALYSIS_STEP_NAME = "02_Azimuth_Mapping"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

PATH_01QC_OUTPUT_NOT_FILTERED = file.path( PATH_EXPERIMENT_OUTPUT, "01_QC", "List_Seurat_Objects_Post_01QC_NotFiltered.rds" )
LITERAL_TITLE = "NK_iPSC_Derived"


MIN.CELLS=3
MIN.FEATURES=200
RIBO_THRESHOLD=10
MITO_THRESHOLD=0
N_GENES_VARIABLES=2000
DO_SCALE=FALSE
DO_CENTER=TRUE
DIM_PCA=50
DIM_UMAP=50
RESOLUTION=0.6
FDRCUTOFF=0.05
FIND_ALL_MARKERS_LOGFC=0.25


PREVENT_OVER_WRITE_REPORT = FALSE

my_palette<-c("#CE4C4B","#343DF4","#E9933D","#AC8B52","#34B7F4")

NUMBER_TOP_SCORING = 20
UMAP_PT_SIZE = 0.8
VISUAL_TRESHOLD_MITO = 10 #Treshold for pct mitochondrial visualisation
VISUAL_TRESHOLD_nFEATURE = 200
VISUAL_TRESHOLD_nCOUNT = 500
LABEL_SIZE = 6 


