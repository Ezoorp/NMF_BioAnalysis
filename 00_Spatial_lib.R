# TODO: FINISH ME!!
# coord_to_df_format(pt_coord) {
#   center <- as.integer(strsplit(center, "x")[[1]])
#   x <- center[1]
#   y <- center[2]
# }

# Convert PT dataframe (two column vectors - 2D coordimates) to Coord Format (XxY) 
pt_df_to_coord_format <- function(pt_df) {
  return(paste0(pt_df[1],"x",pt_df[2]))
}

# remove spots with 0 neighbors from entropy list:
remove_zeroneighbor_spots <- function(spots_coord_list) {
  neighbors_list <- lapply(spots_coord_list, get_1ring, malCells = spots_coord_list)
  neighbors_list <- lapply(neighbors_list, sort)
  
  # Remove 0 neighbor spots:
  zero_neighbor_spots <- which(lapply(neighbors_list, length) == 0)
  if (is_empty(zero_neighbor_spots)) return(spots_coord_list)
  else return(spots_coord_list[-zero_neighbor_spots])
}


# Generate Hex Grid - Only Rings
hex_grid_no_y <- function(hops) {
  # Function to generate points on a hexagonal grid ring with a given distance, 
  # excluding points on the y-axis.
  
  # Create a matrix to store the coordinates of all possible points.
  max_coord <- hops
  coords <- expand.grid(
    x = -max_coord:max_coord, 
    y = -(max_coord * 2) : (max_coord * 2)
  )
  
  # Filter out points based on hexagonal grid rules
  coords <- coords[
    (coords$x %% 2 == 0 & coords$y %% 2 == 0) | 
      (coords$x %% 2 != 0 & coords$y %% 2 != 0) 
    & coords$y != 0, 
  ]
  
  return(coords)
}

# Generate Spots in ring n around a Spot
calculate_hamming_distance <- function(x,y) {
  # Calculate the Hamming distance of a point from the origin (0, 0)
  abs(x) + abs(y)
}
generate_spots <- function(center, distance) {
  # Function to return all points within a given Hamming distance.
  all_points <- hex_grid_no_y(distance)
  all_points$hops <- do.call(calculate_hamming_distance, all_points)
  
  points_at_n_hops <- all_points[all_points$hops == 2*distance, ]
  
  center <- as.integer(strsplit(center, "x")[[1]])
  points_at_n_hops$x <- points_at_n_hops$x + center[1]
  points_at_n_hops$y <- points_at_n_hops$y + center[2]
  
  return(paste0(points_at_n_hops$x,"x",points_at_n_hops$y))   #(points_at_n_hops)
}

# Get Spots from malCells in Ring n
# match(testlist, malCells)

# Randomly Sample 1 spot from Ring n
get_spot_in_ring <- function(center, hops, malCells) {
  mal_ring <- generate_spots(center, hops)
  mal_ring <- malCells[match(mal_ring, malCells)]
  mal_ring <- mal_ring[!is.na(mal_ring)] # remove na's
  
  # if no spots in ring, return NA
  if (length(mal_ring) == 0) return(NA)
  # else return a random spot
  else {
    random_int <- sample(1:length(mal_ring), 1)
    return(mal_ring[random_int])
  }
}

#find spots in 1 hop ring:
# TAKES in a center spot and a list of all cells available in the sample
get_1ring <- function(center, malCells) {
  mal_ring <- generate_spots(center, distance = 1) # 1 hop ring
  mal_ring <- match(mal_ring, malCells)
  mal_ring <- mal_ring[!is.na(mal_ring)] # remove na's
  
  # if no spots in ring, return NA
  if (length(mal_ring) == 0) return(NA)
  else return(mal_ring)
}


# Generate Square Grid 
square_grid <- function(malCells) {
  # Function to generate points on a square grid where each element is the index,
  # index correspond to x/y coord in malCells
  
  # Convert the list into a matrix of coordinates
  coords <- do.call(rbind, strsplit(malCells, "x"))
  coords <- as.data.frame(apply(coords, 2, as.integer))
  colnames(coords) <- c("x", "y")
  
  x_range <- range(coords$x)
  y_range <- range(coords$y)
  grid <- array(NA, dim = c(diff(x_range) + 1, diff(y_range) + 1))
  
  # Fill the grid with indices of coordinates
  for (i in 1:nrow(coords)) {
    grid[coords$x[i] - x_range[1] + 1, coords$y[i] - y_range[1] + 1] <- i
  }
  
  return(list(grid = grid, x_range = x_range, y_range = y_range)) 
}

# Function to check for neighbors in sq grid
check_neighbors_hops <- function(center, grid, x_range, y_range, hops) {
  center <- as.integer(strsplit(center, "x")[[1]])
  x <- center[1]
  y <- center[2]
  x_idx <- x - x_range[1] + 1
  y_idx <- y - y_range[1] + 1
  
  # Generate neighbors at Manhattan distance `hops`
  neighbors <- expand.grid(
    dx = -hops:hops,
    dy = -hops:hops
  )
  
  # Keep only points that are exactly `hops` distance away
  neighbors <- neighbors[abs(neighbors$dx) + abs(neighbors$dy) == hops, ]
  
  results <- sapply(1:nrow(neighbors), function(i) {
    nx <- x_idx + neighbors$dx[i]
    ny <- y_idx + neighbors$dy[i]
    
    if (nx > 0 && ny > 0 && nx <= nrow(grid) && ny <= ncol(grid)) {
      return(grid[nx, ny]) # Returns index or NA
    }
    return(NA)
  })
  
  return(results)
}


# Randomly Sample 1 spot from Square Ring n
get_spot_in_square_ring <- function(center, hops, malCells) {
  grid_params <- square_grid(malCells)
  
  mal_ring <- check_neighbors_hops(center, grid_params$grid, grid_params$x_range, grid_params$y_range, hops)
  mal_ring <- mal_ring[!is.na(mal_ring)] # remove na's
  
  # if no spots in ring, return NA
  if (length(mal_ring) == 0) return(NA)
  # else return a random spot
  else {
    random_int <- sample(1:length(mal_ring), 1)
    return(malCells[mal_ring[random_int]])
  }
}


# Spatial Similarity Analysis over 1 Sample/Seurat Obj
sample_analysis_function <- function(tumor_sample, malProp, spot_dist) {
  cormtx <- matrix(NA, nrow = n_sample_per_ring, ncol = length(test_hops))
  
  colnames(tumor_sample) <- fixSpotNames(colnames(tumor_sample))
  
  # Identify Malig Cells
  tumor_sample$malProp <- malProp
  subset_error <- try(subset(tumor_sample, subset = malProp >= 0.7))
  # ERROR NOT HANDLED CORRECTLY
  if (inherits(subset_error, "try-error")) {
    return(NA)
  } 
  tumor_sample <- subset(tumor_sample, subset = malProp >= 0.7)
  tumor_cell_names <- colnames(tumor_sample)
  n_malig_spots <- ncol(tumor_sample)
  
  print(paste0("Number of Mal Cells in Sample: ", n_malig_spots))
  
  # Skip if too few Malignant Spots
  if ((n_malig_spots < 600) & (spot_dist == 100))  return(NA)
  if ((n_malig_spots < 45) & (spot_dist == 200))  return(NA)
  
  # Standard Preprocessing
  
  tumor_sample <- NormalizeData(tumor_sample, verbose = FALSE)
  tumor_sample <- FindVariableFeatures(tumor_sample, selection.method = "vst", nfeatures = nFeatures, verbose = FALSE)
  malmtx <- GetAssayData(tumor_sample)[VariableFeatures(tumor_sample), ]
  
  # Random Sampling to Find Cor Across Hops
  for (j in seq_along(test_hops)) {
    for (k in 1:n_sample_per_ring) {
      randSpot <- sample(tumor_cell_names, 1)
      
      # 100 vs 200um:
      # TODO: replace with spot dist from tumor sample data structure
      if (spot_dist == 100) {
        spot_in_onehop <- get_spot_in_ring(randSpot, test_hops[j], tumor_cell_names)
      } else {
        spot_in_onehop <- get_spot_in_square_ring(randSpot, test_hops[j], tumor_cell_names)
      }
      if (!is.na(spot_in_onehop)) { # if we find a spot
        if (length(malmtx[,spot_in_onehop]) != nFeatures) print("Found Bug")
        cormtx[k,j] <- cor(malmtx[,randSpot], malmtx[,spot_in_onehop])
      }
    }
  }
  
  
  # Analysis of Correlations:
  n <- colSums(!is.na(cormtx))
  m <- colMeans(cormtx, na.rm = TRUE)
  v <- apply(cormtx, 2, var, na.rm = TRUE)
  
  print(cat("Number of Samples in Each Hop Group:", n))
  return(list(n,m,v))
}



# Spatial Normalized Diversity Analysis over 1 Sample/Seurat Obj
sample_normdiv_analysis_function <- function(malmtx, cormtx) {
  medcor <- med_cor(malmtx)
  # Norm Cor Score
  comp_cormtx <- cormtx / medcor
  # Norm Div Score
  norm_comp_cormtx <- 1 - comp_cormtx
  
  return(norm_comp_cormtx)
}


# Calculate Corr Mtx:
corr_mtx <- function(malmtx) {
  testcor <- cor(as.matrix(malmtx))
  
  # Create Upper Tri Matirx
  testcor[lower.tri(testcor)] <- NA
  # Should appear as upper triangular matrix
  
  # Remove Diagnol (All 1's)
  testcor[lower.tri(testcor,diag=TRUE)] <- NA
  
  return(testcor)
}

# Calculate Median Cor:
median_corr <- function(malmtx) {
  testcor <- corr_mtx(malmtx)
  
  # Median of cor:
  medcor <- median(testcor, na.rm = TRUE)
  # Median because cor is usually skewed
  
  return(medcor)
}


