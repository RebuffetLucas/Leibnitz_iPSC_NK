###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis

ANALYSIS_STEP_NAME = "00_Create_List_Objects"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "iPSC_NK"


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













