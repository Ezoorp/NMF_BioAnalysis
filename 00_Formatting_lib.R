# Making Seurat Obj from Counts Data Files - Visium
#' @param Folder path containing TSV files of counts
#' @return A single seurat object containing the data from each counts file
#' @note the raw counts are usually NOT in tsv file. 
# We have already converted them in a previous script.
load_and_make_seurat_obj <- function(folder) {
  # TEST ME NEED ME: setwd("/Users/levyez/Documents/Professional_Code/")???
  
  # Format into Seurat Obj:
  dataset = list.files(folder);

  if ((!dir.exists(folder)) || (length(dataset[1]) == 0)) {
    stop('Could not find input data')
  }

  # Refactored Code:
  # 1. Load in the counts data from the directory and create Seurat objects
  seurat_objects_list <- lapply(seq_along(dataset), function(i) {
    file_path <- paste0(folder, dataset[i])
    tumor_umi_counts <- read.table(file_path, header = TRUE, row.names = 1, sep = '\t')
    seurat_obj <- CreateSeuratObject(counts = tumor_umi_counts)
    seurat_obj$orig.ident <- NULL
    # Save File Name in Project Name and Sample number in Sample_source
    seurat_obj@project.name <- gsub("*_counts.tsv", "", dataset[i])
    seurat_obj$sample_source <- as.character(i)
    # Modify spot names *before* merging
    colnames(seurat_obj) <- paste0(colnames(seurat_obj),"_",i)

    message(sprintf("obj %d loaded", i))
    return(seurat_obj)
  })

  # 2. Merge Seurat objects
  seurat_object_dataset <- Reduce(function(x, y) JoinLayers(merge(x, y)), seurat_objects_list)
  return(seurat_object_dataset)
} 

#' @param seurat_object_dataset Seurat object containing unfiltered, merged dataset
#' @param seurat_object_dataset_filtered Seurat object containing filtered, merged dataset
#' @param cutoff Percentage of spots remaining after filtering e.g. 50
filter_and_QC_metrics <- function(seurat_object_dataset, seurat_object_dataset_filtered, folder, cutoff = 50) {
  dataset = list.files(folder);
  
  # Calculate spots per sample before and after QC
  preQC_spots_per_sample <- table(seurat_object_dataset$sample_source)
  postQC_spots_per_sample <- table(seurat_object_dataset_filtered$sample_source)

  # Create a data frame for easier manipulation
  spots_per_sample_df <- data.frame(
    Sample = names(preQC_spots_per_sample),
    preQC_SpotCount = as.numeric(preQC_spots_per_sample),
    postQC_SpotCount = as.numeric(postQC_spots_per_sample[names(preQC_spots_per_sample)]) # Match sample names
  )
  spots_per_sample_df$postQC_SpotCount[is.na(spots_per_sample_df$postQC_SpotCount)] <- 0

  # Calculate percentage remaining
  spots_per_sample_df$PercentageRemaining <- (spots_per_sample_df$postQC_SpotCount / spots_per_sample_df$preQC_SpotCount) * 100

  # Identify samples to keep
  samples_to_keep <- spots_per_sample_df$Sample[spots_per_sample_df$PercentageRemaining >= cutoff]

  # Filter the filtered Seurat object
  seurat_object_dataset_filtered <- seurat_object_dataset_filtered[, seurat_object_dataset_filtered$sample_source %in% samples_to_keep]

  # Print QC metrics
  print("QC Metrics:")
  print(spots_per_sample_df)

  # Calculate and print summary metrics
  total_preQC_spots <- sum(spots_per_sample_df$preQC_SpotCount)
  total_postQC_spots <- sum(spots_per_sample_df$postQC_SpotCount[spots_per_sample_df$Sample %in% samples_to_keep])
  total_percentage_remaining <- (total_postQC_spots / total_preQC_spots) * 100
    message(sprintf("Sums (preQC, postQC, %% remaining): %.0f, %.0f, %.1f%%",
          total_preQC_spots, total_postQC_spots,total_percentage_remaining))

  # Split into list of Seurat objects
  sample_names <- unique(seurat_object_dataset_filtered$sample_source)
  tumor_seurat_objs <- lapply(sample_names, function(sample_name) {
    message(sprintf('Create QC\'d seurat_obj for sample %s ...', sample_name))
    subset(seurat_object_dataset_filtered, subset = sample_source == sample_name)
  })
  names(tumor_seurat_objs) <- sample_names

  # Attach sample names to each sample's seurat object in @ project name
  for (i in seq_along(tumor_seurat_objs)) {
    dataset_index <- as.integer(names(tumor_seurat_objs)[i])
    if (!is.na(dataset_index) && dataset_index <= length(dataset)) {
      tumor_seurat_objs[[i]]@project.name <- paste0(gsub("*_counts.tsv", "", dataset[dataset_index]))
    }
  }

  return(tumor_seurat_objs)
}
