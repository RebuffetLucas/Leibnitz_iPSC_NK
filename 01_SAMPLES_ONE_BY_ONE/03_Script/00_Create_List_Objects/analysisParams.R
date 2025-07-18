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
SCALE_FACTOR = 10000
DO_SCALE=TRUE
DO_CENTER=TRUE
DIM_PCA=50
DIM_UMAP=50
RESOLUTION=0.6
FDRCUTOFF=0.05
FIND_ALL_MARKERS_LOGFC=0.25

SEED = 19

PREVENT_OVER_WRITE_REPORT = FALSE
NUMBER_TOP_SCORING = 20


#Filter DEG
MARKERS_PVALUE_TRESHOLD = 5e-2


#Graphical Params
PT_SIZE_FEATURE_SCATTER = 0.5
PT_SIZE_FEATURE_PLOT = 1.2

my_palette<-c("#CE4C4B","#343DF4","#E9933D","#AC8B52","#34B7F4")

MAX_CUTOFF_SCORE_MODULE = 1.5
LABEL_SIZE = 6 
UMAP_PT_SIZE = 0.8
NUMBER_PCs_CORRELATION_PLOTS = 30


#Visual tresholds for QC filtering dataviz'
VISUAL_TRESHOLD_MITO = 10 #Treshold for pct mitochondrial visualisation
VISUAL_TRESHOLD_nFEATURE = 200
VISUAL_TRESHOLD_nCOUNT = 500



#List Genes of interest for ploting
LIST_GENE_1 = c("FCER1G", "SPON2", "GZMB",
                "XCL1", "XCL2", "GZMK",
                "CD3E", "CD3D", "KLRC2")

LIST_GENE_2 = c("NKG7", "LYZ", "PRF1",
                "CD44", "IL7R", "CD2",
                "IL32", "CCL5", "PRDM1")




