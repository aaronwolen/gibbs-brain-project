#' Calculate correlations for paired features across expressionSet objects 
#'   
#' @param x expressionSet object
#' @param y expressionSet object
#' @param feature.key a data.frame or matrix with three columns: the first
#'   contains a grouping variable, the second and third contain feature
#'   indentifiers in eSet x and y, respectively. Correlations are calculated for
#'   all combinations of x and y features, within each group.
#' @param p.adjust.method method of p-value adjustment
#' @param ... additional arguments passed to cor
#' 
#' @return
#' A data.frame containing the values from feature.key as well as the 
#' correlation statistic ('r'), the number of non-NA cases used in the 
#' calculation ('n'), the p-value of the test ('pvalue') and adjusted p-value of
#' the test ('adj.pvalue').


cor_cross_eset <- function(x, y, feature.key, p.adjust.method = "BH", ...) {
  require(Biobase, quietly = TRUE)
  
  p.adjust.method <- match.arg(p.adjust.method, p.adjust.methods)
  
  # TODO: Check feature.key
  
  # Capture feature.key column names
  labels <- names(feature.key)
  names(feature.key) <- c("group", "x", "y")
  
  # Check eSet objects
  if( any(c(class(x), class(y)) != "ExpressionSet") ) {
    stop("Both x and y must be eSet objects.", call. = F)
  }
  
  if( ncol(x) != ncol(y) ) {
    stop("Both eSets must contain the same number of samples.", call. = F)
  }
  
  # Check samples
  samples <- intersect(sampleNames(x), sampleNames(y))
  
  if( length(samples) != ncol(x) ) {
    stop("Both eSets must contain identical samples.", call. = F)
  }
  
  # Ensure samples are aligned
  x <- x[, samples]
  y <- y[, samples]
  
  # Extract transposed expression matrices
  x.mat <- t(exprs(x))
  y.mat <- t(exprs(y))
  
  # Results data.frame
  cors.df <- list()
  
  # Loop through each level of grouping variable
  for(g in unique(feature.key$group)) {
    
    # Probes from each eSet within current group
    xp <- unique(feature.key$x[feature.key$group == g])
    yp <- unique(feature.key$y[feature.key$group == g])
    
    out <- cor(x.mat[, xp], y.mat[, yp], use = "pairwise.complete", ...)
    
    if(length(xp) == 1 | length(yp) == 1) {
      out <- data.frame(Var1 = xp, Var2 = yp, Freq = as.numeric(out))
    } else {
      out <- as.data.frame.table(out)  
    }
    
    # Count number of values used in each correlation
    out$n <- apply(as.matrix(out), 1, function(x)
      sum(!is.na(x.mat[, x[1]]) & !is.na(y.mat[, x[2]])))
    
    cors.df[[g]] <- data.frame(group = g, out)
  }
  
  # Unlist
  cors.df <- data.frame(do.call("rbind", cors.df), row.names = NULL)
  
  # Apply labels captured from feature.key
  cors.df <- cors.df[, c(1:3, 5, 4)]
  names(cors.df) <- c(labels, "n", "r")
  
  # Calculate p-values (using code adapted from WGCNA)
  cors.df$pvalue <- with(cors.df, 
    2 * pt(abs(sqrt(n - 2) * r/sqrt(1 - r^2)), n - 2, lower.tail = FALSE))
  
  # Calculate adjusted pvalues
  cors.df$adj.pvalue <- p.adjust(cors.df$pvalue, method = p.adjust.method)
  
  return(cors.df)
}

