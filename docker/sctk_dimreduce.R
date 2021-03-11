suppressMessages(suppressWarnings(library(singleCellTK)))

# Importing CellRanger data

# Auto determines V2 vs. V3
# SCTK expects a specific folder structure
# /data_directory/cellranger_sample_folder/outs/filtered_feature_bc_matrix
# Insided filtered_feature_bc_matrix (gzipped or not)
#     barcodes.tsv
#     features.tsv
#     matrix.mtx
# ex: ./data/5k_pbmc/outs/filtered_feature_bc_matrix
#sce <- importCellRanger(
#    cellRangerDirs = "./data",
#    sampleDirs = "5k_pbmc",
#    sampleNames = "5k_pbmc",
#    dataType = "filtered",
#    gzipped = "auto"
#)


# args from command line:
args <- commandArgs(TRUE)
RAW_COUNT_MATRIX <- args[1]
OUTPUT_UMAP_BASE <- 'umap_matrix'

# Import counts as a 
counts <- read.table(
    file = RAW_COUNT_MATRIX,
    sep = "\t",
    row.names = 1
)

# Create an SCE object from the counts
sce <- SingleCellExperiment(
    assays=list(counts=counts)
)

# PCA as a pre-processing step prior to UMAP dimensionality reduction.
# Important to decorrelate the counts prior to UMAP.
# Also log normalizes the counts.
# Default to 50 PCA dims. Higher than typically needed.
# Performing UMAP (on PCA) and adding UMAP data to the SCE object
sce <- getUMAP(
    inSCE = sce, 
    useAssay = "counts", 
    reducedDimName = "UMAP",
    logNorm = TRUE,
    nNeighbors = 30, 
    nIterations = 200, 
    alpha = 1,
    minDist = 0.01,
    pca = TRUE,
    initialDims = 50
)

# Export UMAP values to file
df.umap <- SingleCellExperiment::reducedDim(sce, "UMAP")