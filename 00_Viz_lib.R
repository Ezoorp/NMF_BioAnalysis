# TODO: make spatial plots more scrunched up (smaller y height)
# generally looks better / more clear pattern changes with less white space btw pts

# Input: Vector with Coordinates as Names


# Ultrasimple 1 Variable plotter
# Input: vector with coordinates in name (format: X#xY# or XX#xY#)
# Set imagePath = NA if printing to IDE
fastSpatialPlot <- function(visualVector, 
                            imagePath = NULL, 
                            DiscreteFlag = FALSE, 
                            Groups = NULL, 
                            RotateFlag = FALSE,
                            Title = NULL,
                            color_high = "red",
                            color_mid = "white",
                            color_low = "blue",
                            legend_title = "") {
  
  # Check if named coordinates, if so fix them:
  if(all(grepl("^X\\d+x\\d+$", names(visualVector)))) {
    names(visualVector) <- fixSpotNames(names(visualVector))
  }
  
  coord <- t(matrix(as.numeric(unlist(strsplit(names(visualVector),"x"))),nrow=2))
  x = coord[,1]
  y = coord[,2]
  
  if (RotateFlag == TRUE){
    x = coord[,2]
    y = max(coord[,1]) - coord[,1]
  }
  
  fig.df <- data.frame(
    x=x,
    y=y,
    value=visualVector
  )
  # TODO: change to ggsave with large save?
  if(!is.null(imagePath)) png(filename = imagePath, width=700, height=500)
  if (DiscreteFlag == TRUE) {
    p <- ggplot(fig.df,aes(x=x,y=y, color = as.factor(value))) +
      geom_point() + 
      scale_color_manual(values = DiscretePalette(Groups)) + 
      ggtitle(Title)
    print(p)
  } else {
    p <- ggplot(fig.df,aes(x=x,y=y, color = value)) +
      geom_point() +
      scale_color_gradient2(high = color_high, mid = color_mid, low = color_low, midpoint = mean(visualVector)) +
      ggtitle(Title) +
      labs(colour = legend_title)
      #theme(panel.background = element_rect(fill = 'black', colour = 'white'))
    print(p)
  }
  if(!is.null(imagePath)) dev.off()
}


# Fixing rownames of post QC data from X#1x#2_#3 to #1x#2 
fixSpotNames <- function(rownames) {
  rownames <- gsub("X|_.*", "", rownames)
  rownames
} 


# Metaprogram Plotting
# Keeps each Metaprogram the same color for consistency accross plots.
MPfastSpatialPlot <- function(visualVector, 
                              imagePath = NULL, 
                              RotateFlag = FALSE,
                              Title = NULL,
                              Height = 5,
                              Width = 11) {
  
  coord <- t(matrix(as.numeric(unlist(strsplit(names(visualVector),"x"))),nrow=2))
  x = coord[,1]
  y = coord[,2]
  
  if (RotateFlag == TRUE){
    x = coord[,2]
    y = max(coord[,1]) - coord[,1]
  }
  
  fig.df <- data.frame(
    x=x,
    y=y,
    value=visualVector
  )
  
  # Make MP Color Palette
  MP_color_palette <- DiscretePalette(36)
  MP_color_palette[2] <- "#FFFF00"
  MP_color_palette[37:38] <- c("#AEB8B7","#ED7480")
  names(MP_color_palette) <- paste0("MP", 1:38)
  # Get the unique MP values present in the data
  unique_mps <- unique(visualVector)
  color_mapping <- MP_color_palette[unique_mps]
  
  # TODO: change to ggsave with large save?
  
  
  p <- ggplot(fig.df,aes(x=x,y=y, color = as.factor(value))) +
    geom_point() + 
    scale_color_manual(values = color_mapping) + 
    ggtitle(Title) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,1),"in"))
  
  if(!is.null(imagePath)) ggsave(filename = imagePath, 
                                 plot = p, 
                                 units = "in", 
                                 height = Height,
                                 width = Width)
}
