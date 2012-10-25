#' Cross expressionSet scatterplots
#' 
#' @param x an ExpressionSet object
#' @param y an ExpressionSet object
#' @param feature.x character string identifying a feature in x
#' @param feature.y character string identifying a feature in y
#' @param mapping aesthetics passed for each layer with aes()
#' @param id.name character string identifying a phenoData variable that can be
#'   used to identify common samples in both ExpressionSet ojects
#' @param facets faceting formula to use. (see ?qplot for more information)

plot_cross_eset <- function(x, y, feature.x, feature.y, mapping = aes(), 
  id.name = NULL, facets = NULL) {
  
  suppressPackageStartupMessages(require(ggplot2, quietly = TRUE))
  require(Biobase, quietly = TRUE)
  require(ggaffy, quietly = TRUE)
  
  stopifnot(all(c(class(x), class(y)) %in% "ExpressionSet"))
  
  # Can only handle one feature from each eSet object
  stopifnot(length(c(feature.x, feature.y)) == 2)
  
  # Assemble vector of variables
  if(missing(id.name)) id.name <- "sample"
  aes.vars <- unique(as.character(unlist(mapping)))
  facet.vars <- all.vars(facets)
  vars <- c("sample", "expression", id.name, aes.vars, facet.vars)
  
  # Extract and combine eSet data.frames
  exprs.df <- rbind(
    get_eset_data(x, feature.x, shape = "long")[, vars],
    get_eset_data(y, feature.y, shape = "long")[, vars]
  )
  exprs.df$eset <- c(rep("x", ncol(x)), rep("y", ncol(y)))
  
  # Verify data points are uniquely identified
  if (any(table(exprs.df[, c("eset", id.name, facet.vars)]) > 1)) {
    stop(id.name, "does not uniquely identify samples.", call. = FALSE)
  }
     
  # Format data for plotting
  cast.f <- as.formula(paste(
    paste(c(id.name, aes.vars, facet.vars), collapse = "+"), 
    "eset", sep = "~")
  )
  exprs.df <- reshape2::dcast(exprs.df, cast.f, value.var = "expression")
  
  p <- ggplot(exprs.df, aes(x = x, y = y)) + geom_point(mapping)
  
  if (is.null(facets)) {
    p <- p + facet_null()
  } else if (is.formula(facets) && length(facets) == 2) {
    p <- p + facet_wrap(facets)
  } else {
    p <- p + facet_grid(facets = deparse(facets))
  }
  
  p + xlab(feature.x) + ylab(feature.y)
}
