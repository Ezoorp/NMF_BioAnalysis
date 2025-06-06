# Load Libraries:
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(Rcpp))
setwd("/Users/levyez/Documents/Professional_Code/NMF_BioAnalysis/")

source('00_Viz_lib.R')
source('00_Formatting_lib.R')
source('00_NMF_lib.R')
sourceCpp('00_NMF_fast_lib.cpp')

output_dir = "Outputs/"
dir.create(output_dir, showWarnings = "FALSE")

##########################################
## Load in Raw Counts Data from directory:
# By hand on laptop:
#folder <- "Wu_2022_Liver"
# By hand on HPC:
folder <- "/data/CDSL_MAlab/levyez/00_cell_states/Visium_test_data/Wu_2022_Liver"
# With Bash + HPC Script:
#args <- commandArgs(trailingOnly = TRUE)
#folder=args[1]

path.folder <- paste0(folder,"/")
print(paste("Input folder:", folder))

seurat_object_dataset <- load_and_make_seurat_obj(path.folder)

## Quality Control:
# Filer out low count gene, features, and high Mito spots
seurat_object_dataset[["percent.mt"]] <- PercentageFeatureSet(seurat_object_dataset, pattern = "^MT-")

seurat_object_dataset_filtered <- subset(seurat_object_dataset, subset =
                                   nCount_RNA > 1000 &
                                   nFeature_RNA > 400 &
                                   percent.mt < 40)

# Remove Low Quality Samples, Reformat into List of Seurat Objects
seurat_objs_list <- filter_and_QC_metrics(seurat_object_dataset, seurat_object_dataset_filtered, path.folder)

# Save Data in RDS File
saveRDS(seurat_objs_list, file = paste0(output_dir,folder,'_PostQC_Seurat_Objects.RData'))
message(sprintf("Seurat Objects saved to path: %s", paste0(output_dir,folder,'_PostQC.RData')))


##########################################
## Normalization + Highly Variable Features:
# TO Normalize:
# 1. Divide each cell by the total number of molecules measured in the cell
# 2. Multiply that number by a scaling factor (i.e. 10000)
# 3. Add 1, and take a natural log

hvg_features <- 2000

seurat_objs_list <- lapply(seurat_objs_list, function(obj){ 
  message(sprintf('Normalizing + HVG sample %s, of total %d ',
                  obj$sample_source[[1]], length(seurat_objs_list)))
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = hvg_features)
  obj
})

##########################################
## Begin SpaCET:

##########################################
## Begin inferCNV w. SpaCET identified Normal:

##########################################
## Nonnegative Matrix Factorization:
#!!!!!!! USER SETTINGS:
# 1.Choose whether to use ALL cells or MALIGNANT cells:
malig_bool <- FALSE
# 2. Setting the number of components (K) we want NMF to identify:
rank_test <- 2:4 # K equal to 2 to 11
# 3. Set # of genes to represent each component
max_program_genes <- 50

# File Naming:
nmf_name <- sprintf('k%d_%d.hvg%d', 
                    min(rank_test), max(rank_test), hvg_features, malig_bool)

path.nmf.folder <- paste0(output_dir, "NMF_",nmf_name)
if (!dir.exists(path.nmf.folder))
  dir.create(path.nmf.folder)

path.nmf.file <- file.path(path.nmf.folder, paste0('NMF_raw_output_', nmf_name, '.rds'))
###############
# RUN NMF ALGORITHM:
if (!file.exists(path.nmf.file)){
  geneNMF.programs <- multiSampleNMF(
    seurat_objs_list, 
    k = rank_test, 
    nfeatures = hvg_features)
  
  saveRDS(geneNMF.programs, file = path.nmf.file)
} else {
  geneNMF.programs <- readRDS(file = path.nmf.file)
}


#View(geneNMF.programs)

###############
# Plot Reconstruction Error from NMF
plot_nmf_tolerances(geneNMF.programs, path.nmf.folder, nmf_name)

###############
# Identify Top Weight Genes in Each Component
message("Max # of Gene per Programs: ", max_program_genes)
nmf.genes <- getNMFgenes(nmf.res=geneNMF.programs, max.genes=max_program_genes)

# SAVE NMF Components
path.nmf.genes.file <- paste0(path.nmf.folder, '/NMF_components_', nmf_name, '.rds')
saveRDS(nmf.genes, file = path.nmf.genes.file)

###############
## Filter NMF components (may be called programs in some cases)
# Set cutoff for intrasample similarity, intrasample redundancy, and intersample similarity:
intrasample_similarity_cutoff <- 0.6
intrasample_redundancy_cutoff <- 0.2
intersample_similarity_cutoff <- 0.2
nprogs <- length(nmf.genes)
message("Total # of initial Programs: ", nprogs)

# calculate the similarities within each samples and do filtering
min_intra_sim_robust <- max_program_genes * intrasample_similarity_cutoff
max_intra_sim_redundant <- max_program_genes * intrasample_redundancy_cutoff
min_inter_sim_robust <- max_program_genes * intersample_similarity_cutoff

# Create J matrix (Intersection Matrix)
J <- calculateOverlap(nmf.genes)

# Heatmap of Initial J matrix:
path.intersectionmatrix.heatmap <- paste0(path.nmf.folder, '/NMF_interesection_',nmf_name,'.png')
plotIntersectionMatrix(J, output_file = TRUE, path.intersectionmatrix.heatmap)

# 1. filter programs that are similar intrasamplewise
robust.intra.progs <- filterIntraSimilarProgs(
  Jmatrix = J,
  nmf.genes = nmf.genes, 
  ranks = rank_test, 
  sample_names = names(seurat_objs_list), 
  min_intra_sim_robust = min_intra_sim_robust
)

# DEBUG: Programs Per Sample:
# table(gsub('.k\\d+\\.p\\d+', '', robust.intra.progs))


# SAVE Filtered NMF Components
path.nmf.filtered.genes.file <- paste0(path.nmf.folder, '/Filtered_NMF_components_', nmf_name, '.rds')
saveRDS(robust.intra.progs, file = path.nmf.filtered.genes.file)

############################
# NEXT STEPS: Run the Inter Sample Similarity + Clustering accross all dataset

## OR: Run Inter Sample Similarity for THIS dataset -> clustering + analysis
############################

# 2. filtering across samples to remove non-robust program 
# Constrain J matrix to just remaining programs:
J_filtered <- J[robust.intra.progs,robust.intra.progs]

# TODO: rewrite code with lapply / avoid For loops
robust.inter.progs.intersection <- filterInterSimilarProgs(
  Jmatrix = J_filtered,
  robust.intra.progs = robust.intra.progs, 
  min_inter_sim_robust = min_inter_sim_robust
)
robust.inter.progs <- names(robust.inter.progs.intersection)

# 3. remove redundant programs within each sample
# Constrain J matrix to just remaining programs:
J_filtered <- J_filtered[robust.inter.progs,robust.inter.progs]

redunt_progs <- c()
keep_progs <- c()

for (i in seq_along(robust.inter.progs.intersection)){
  cur_prog <- names(robust.inter.progs.intersection)[i]
  if (cur_prog %in% redunt_progs) 
    next
  keep_progs <- c(keep_progs, cur_prog)
  sample_i <- gsub('.k\\d+\\.p\\d+', '', cur_prog)
  if (i >= length(robust.inter.progs.intersection)) next
  for (j in (i + 1) : length(robust.inter.progs.intersection)){
    tmp_prog <- names(robust.inter.progs.intersection)[j]
    if (tmp_prog %in% redunt_progs) next
    sample_j <- gsub('.k\\d+\\.p\\d+', '', tmp_prog)
    if (sample_i == sample_j){
      if (J[cur_prog, tmp_prog] > max_intra_sim_redundant)
        redunt_progs <- c(redunt_progs, tmp_prog)
    }
  }
}



###########################
# Vizualizing Programs after ALL Filtering:
# Constrain J matrix to just remaining programs:
J_filtered <- J_filtered[keep_progs,keep_progs]

path.filtered.intersectionmatrix.heatmap <- paste0(path.nmf.folder, '/NMF_filtered_interesection_',nmf_name,'.png')
plotIntersectionMatrix(J_filtered, output_file = TRUE, path.filtered.intersectionmatrix.heatmap)

