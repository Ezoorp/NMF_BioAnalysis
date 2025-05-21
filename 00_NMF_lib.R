


#' Run NMF on a list of Seurat objects
#'
#' Given a list of Seurat objects, run non-negative matrix factorization on 
#' each sample individually, over a range of target NMF components (k).
#'
#' @param obj.list A list of Seurat objects
#' @param k Number of target components for NMF (can be a vector)
#' @param assay Get data matrix from this assay
#' @param slot Get data matrix from this slot (=layer)
#' @param nfeatures Number of HVG, if calculate_hvg=TRUE
#' @param L1 L1 regularization term for NMF
#' @param min.cells.per.sample Minimum numer of cells per sample (smaller 
#'     samples will be ignored)
#' @param min.exp Minimum average log-expression value for retaining genes
#' @param max.exp Maximum average log-expression value for retaining genes
#' @param seed Random seed     
#'     
#' @return Returns a list of NMF programs, one for each sample and for each
#'     value of 'k'. The format of each program in the list follosw the
#'     structure of \code{\link[RcppML]{nmf}} factorization models.
#'
#' @examples
#' library(Seurat)
#' data(sampleObj)
#' geneNMF_programs <- multiNMF(list(sampleObj), k=5)
#' 
#' @importFrom RcppML nmf
multiSampleNMF <- function(obj.list, assay="RNA", slot="data", k=5:6,
                           nfeatures = 2000, L1=c(0,0),
                           min.exp=0.01, max.exp=3.0,
                           min.cells.per.sample = 10,
                           seed = 1234
                           ) {
  
    #set.seed(seed)
  
    nmf.res <- lapply(obj.list, function(obj){
    message(sprintf('Calculating NMF of sample %s, total %d ...', 
                    obj$sample_source[[1]], length(obj.list)))
    mat <- as.matrix(GetAssayData(subset(obj,features = VariableFeatures(obj)), assay="RNA"))

    
    res.k <- lapply(k, function(k.this){
      message(
        sprintf('Calculating NMF of sample %s, K %d, minK %d, maxK %d ...', 
                obj$sample_source[[1]], k.this, min(k), max(k)))
      model <- RcppML::nmf(mat, k = k.this, L1 = L1, verbose = FALSE) #, seed = seed)
      rownames(model$h) <- paste0("pattern", 1:nrow(model$h))
      colnames(model$h) <- colnames(mat)
      rownames(model$w) <- rownames(mat)
      colnames(model$w) <- paste0("pattern",1:ncol(model$w))
      model
    })
    #names(res.k) <- paste0(obj$sample_source[[1]], ".k",k)
    names(res.k) <- paste0("k",k)
    res.k
  })
  nmf.res <- unlist(nmf.res, recursive = FALSE)
  
  return(nmf.res)
}  

#' Plot NMF Tolerance Histogram
#' Given a list of NMF programs, plot a histogram of the tolerance values.
#'
#' @param geneNMF.programs A list of NMF programs (output from multiSampleNMF)
#' @param path.nmf.folder The directory to save the histogram to
#' @param nmf_name Name of the NMF run, used for the filename
#'
#' @return Saves a PNG file of the histogram to the specified directory.
plot_nmf_tolerances <- function(geneNMF.programs, path.nmf.folder, nmf_name) {
  # Extract tolerance values
  tolerances <- sapply(geneNMF.programs, function(x) x$tol)

  # Create the histogram
  hist_path <- paste0(path.nmf.folder, "/NMF_tolerance_histogram_",nmf_name,".png")
  png(hist_path, width = 800, height = 600) # Open PNG device
  hist(tolerances, breaks = 20, main = "Histogram of NMF Tolerances",
       xlab = "Tolerance", col = "skyblue", border = "steelblue")
  dev.off() # Close PNG device
  message(sprintf("NMF tolerance histogram saved to: %s", hist_path))
}


#' Get list of genes for each NMF program
#'
#' Run it over a list of NMF models obtained using \code{multiNMF()}
#'
#' @param nmf.res A list of NMF models obtained using \code{multiNMF()}
#' @param method Parameter passed to \code{\link[NMF]{extractFeatures}} to
#'     obtain top genes for each program. When 'method' is a number between 0 
#'     and 1, it indicates
#'     the minimum relative basis contribution above which the feature is
#'     selected, i.e. how specific is a gene for a given program.
#' @param max.genes Max number of genes for each program 
#'     
#' @return Returns a list of top genes for each gene program found
#'     by \code{multiNMF()}
#' @importFrom NMF extractFeatures
getNMFgenes <- function(nmf.res, method=0.5, max.genes=200) {
  nmf.genes <- lapply(nmf.res, function(model) {
    
    emb <- model$h
    load <- model$w
    
    m <- NMF::extractFeatures(load, method=method)
    m <- lapply(m, function(x){
      genes <- rownames(load)[x]
      head(genes, min(length(genes), max.genes))
    })
    
    names(m) <- paste0("p",seq(1,length(m)))
    m
  })
  
  nmf.genes <- unlist(nmf.genes, recursive = FALSE)
  return(nmf.genes)
}



# ' Calculate J matrix (intersection matrix) for NMF genes
#' @param nmf.genes A list of NMF genes
#' @return A matrix representing the intersection of NMF genes
#' @import RCPP and calculateOverlapOptimized
calculateOverlap <- function(nmf.genes) {
  J <- matrix(data=0, ncol=length(nmf.genes),
              nrow = length(nmf.genes))

  # Running CPP J (Intersection) Matrix
  st.time <- Sys.time()
  J <- calculateOverlapOptimized(nmf.genes)
  message('Calculating Intersection Matrix with CPP:')
  print(Sys.time() - st.time)

  # Label J matrix
  colnames(J) <- names(nmf.genes)
  rownames(J) <- names(nmf.genes)
}

#' Filter NMF programs based on intra-sample similarity
#' @param Jmatrix The J matrix (intersection matrix) calculated from the NMF genes
#' @param nmf.genes The list of NMF genes
#' @param ranks The ranks of the NMF genes
#' @param sample_names The names of the samples
#' @param min_intra_sim_robust The minimum intra-sample similarity threshold
#' @param verbose Whether to print verbose output
#' @return A vector of filtered program names
filterIntraSimilarProgs <- function(
    Jmatrix, nmf.genes, ranks, sample_names, min_intra_sim_robust, verbose=FALSE
){
  filtered_program <- c()
  sample_nprogs <- sum(ranks)
  st.time <- Sys.time()
  for (i in seq_along(sample_names)){
    if (verbose)
      print(sprintf('Calculating intra-sample similarity of sample %s, %d of %d', 
                    sample_names[i], i, length(sample_names)))
    idx.prog.st <- (i - 1) * sample_nprogs + 1
    idx.prog.en <- i * sample_nprogs
    
    sample.progs <- names(nmf.genes)[idx.prog.st: idx.prog.en]
    Jsample <- Jmatrix[sample.progs, sample.progs]
    
    for (pi in 1:sample_nprogs){
      if (verbose)
        print(sprintf('Calculating intra-sample similarity of sample program %s, %d of %d', 
                      sample_names[i], pi, sample_nprogs))
      tmpJ <- Jsample[pi, c(seq(1, sample_nprogs)[seq(1, sample_nprogs) != pi])]
      tmpMax <- max(tmpJ)
      if (tmpMax >= min_intra_sim_robust)
        filtered_program <- c(filtered_program, sample.progs[pi])
    }
  }
  print('Calculating and filter intra-Similarities...')
  print(Sys.time() - st.time)
  message("Total # of intra_similar filtered Programs: ", length(filtered_program))
  filtered_program
}


# Plot Heatmap of Initial J matrix:
plotIntersectionMatrix <- function(J, output_file = TRUE, path.nmf.folder, nmf_name) {
  if (output_file) {
    path.intersectionmatrix.heatmap <- paste0(path.nmf.folder, '/NMF_interesection_',nmf_name,'.png')
    png(path.intersectionmatrix.heatmap, 
      width=8.5,
      height=8,
      units="in",
      res=1200)
  }
  col_fun = colorRamp2(c(0, 25, 50), c("blue", "white", "red"))
  Heatmap(J, col = col_fun,
          heatmap_legend_param = list(
            title = "Intersection", at = c(0, 25, 50), 
            labels = c("0", "25", "50")
          )
          )
  if (output_file) { 
    dev.off()
    message(sprintf("Intersection matrix heatmap saved to: %s", path.intersectionmatrix.heatmap))
  }
}