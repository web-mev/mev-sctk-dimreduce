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
sce <- importCellRanger(
    cellRangerDirs = "./data",
    sampleDirs = "5k_pbmc",
    sampleNames = "5k_pbmc",
    dataType = "filtered",
    gzipped = "auto"
)

# Adding UMAP data to the SCE object
sce <- getUMAP(
    inSCE = sce, 
    useAssay = "counts", 
    reducedDimName = "UMAP"
)

df.umap <- SingleCellExperiment::reducedDim(sce, "UMAP")