#' Heatmap visualization of small matrices
#' 
#' @param mat numerical matrix with dimnames
#' @param tile.labels logical, if TRUE each tile's numerical value is provided
#' @param symmetric logical, if TRUE values along the diagnol are ignored and
#'   used as labels instead of axis text

heat_matrix <- function(mat, tile.labels = TRUE, symmetric = FALSE,
  legend.title = "value", font.size = 5) {
  require(ggplot2, quietly = TRUE)
  
  df <- as.data.frame.table(mat)
  
  df$Var1 <- factor(df$Var1, levels = rownames(mat))
  df$Var2 <- factor(df$Var2, levels = colnames(mat))
  
  p <- ggplot(df) + aes(x = Var1, y = Var2) +
    geom_tile(aes(fill = Freq))
  
  p <- p + scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))
    
  if(symmetric) {
    diag.rows <- which(df$Var1 == df$Var2)
    df$Freq[diag.rows] <- NA
    p <- p %+% df
    p <- p + 
      geom_text(data = df[diag.rows, ], 
        aes(label = Var1), color = "white", size = font.size) +
      opts(axis.text.x = theme_blank(), axis.text.y = theme_blank())
  }
  
  if(tile.labels) {
    p <- p + 
      geom_text(data = df[!is.na(df$Freq),],
        aes(label = round(Freq, d = 2)), color = "white", size = font.size)
  }
  
  p <- p + opts(axis.title.x = theme_blank(), axis.title.y = theme_blank(),
    axis.ticks = theme_blank())
  
  p <- p + guides(fill = guide_colorbar(title = legend.title))
  
  return(p)
}
