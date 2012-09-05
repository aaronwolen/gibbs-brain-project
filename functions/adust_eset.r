#' Adjust expressionSet objects for supplied covariates
#'
#' @param eset expressionSet object
#' @param formula formula to be passed to lm_eset
#' @param output path of the output file, if NULL adjusted eSet object is returned
#' @param replace.existing logical, should existing output path be overwritten? 
#' 
#' Expression data is adjusted independently for each tissue. 

adjust_eset <- function(eset, formula, output = NULL, replace.existing = FALSE) {
  
  require(Biobase, quietly = TRUE)
  source("functions/lm_eset.r")
  
  # Check for existing output file
  if( file.exists(output) & !replace.existing) {
    warning(cat(output, "already exists.\n"), call. = FALSE)
    return(NULL)
  }
  
  # Parse formula
  lm.f <- as.formula(match.call()[["formula"]])
  lm.vars <- all.vars(lm.f)
  
  # Loop through each tissue and perform adjustment
  exp.mat <- list()
  for(t in levels(eset$tissue)) {    
    exp.mat[[t]] <- do.call(lm_eset, list(lm.f, eset = eset[, eset$tissue == t],
      use = "pairwise.complete", resid.only = TRUE))
  }
  
  exp.mat <- do.call("cbind", exp.mat)
  
  stopifnot(all(rownames(exp.mat) == rownames(exprs(eset))))
  stopifnot(all(colnames(exp.mat) == colnames(exprs(eset))))
  
  exprs(eset) <- exp.mat
  
  # Add annotation to eset object
  msg <- paste(lm.vars, collapse = ", ")
  msg <- sub(paste(",", tail(lm.vars, 1)), paste(" and", tail(lm.vars, 1)), msg)
  
  eset@annotation <- c(eset@annotation,
    paste("Adjusted data for covariates", msg))
  
  # Return results
  if(is.null(output)) {
    return(eset)
  }
  
  # Assign adjusted eSet object to provided variable name and export
  eset.name <- as.character(match.call()[["eset"]])
  assign(eset.name, eset)
  save(list = eset.name, file = output)
}
