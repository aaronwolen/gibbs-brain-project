# plotmatrix function that accepts a data.frame with additional variables for
# aes mapping

plot_matrix <- function (data, facet.vars, mapping = aes(), 
  size = 1, alpha = 1, density = FALSE, bins = 30) {

  if(missing(facet.vars)) {
    facet.vars <- colnames(data)
  }
  
  if(length(mapping) > 0) {
    aes.data <- data[, setdiff(names(data), facet.vars)]
  }
  
  data <- data[, facet.vars]
  
  grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
  grid <- subset(grid, x != y)
  
  # data.frame with xy coordinates
  all <- lapply(1:nrow(grid), function(i) {
    xcol <- grid[i, "x"]
    ycol <- grid[i, "y"]
    data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol], 
               x = data[, xcol], y = data[, ycol], data)
  })
  
  # add aes variables
  if(length(mapping) > 0) {
    all <- lapply(all, cbind, aes.data)
  }
  
  all <- do.call("rbind", all)
  
  all$xvar <- factor(all$xvar, levels = names(data))
  all$yvar <- factor(all$yvar, levels = names(data))
  
  densities <- do.call("rbind", lapply(1:ncol(data), function(i) {
    data.frame(xvar = names(data)[i], yvar = names(data)[i], 
               x = data[, i])
  }))
  
  xy.mapping <- aes_string(x = "x", y = "y")
  class(xy.mapping) <- "uneval"
  
  p <- ggplot(all, xy.mapping) + facet_grid(xvar ~ yvar, scales = "free")
  
  if(density) {
    p <- p + stat_binhex(bins = bins, geom = "hex")
  } else {
    if("size" %in% names(mapping)) {
      p <- p + geom_point(mapping, na.rm = TRUE, alpha = alpha)
    } else {
      p <- p + geom_point(mapping, na.rm = TRUE, size = size, alpha = alpha)  
    }
  }
    
  p + stat_density(data = densities, 
      aes(x = x, y = ..scaled.. * diff(range(x)) + min(x)), 
      position = "identity", colour = "grey20", geom = "line") +
    opts(axis.title.x = theme_blank(), axis.title.y = theme_blank())
}
