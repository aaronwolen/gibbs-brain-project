#' Fitting linear models across ExpressionSet objects
#' 
#' @param formula an object of class "formula"
#' @param eset an ExpressionSet object
#' @param probesets character vector of probesets to be used in the fitting
#' @param use default is to use only "complete" samples (ie, no missing values),
#' "pairwise.complete" removes missing values on a feature by feature basis and
#' returns a list of individual lm objects

#' @param ... additional arguments passed to lm

#' This currently returns different kinds of results depending on whether use
#' equals complete or pairwise complete.

lm_eset <- function(formula, eset, probesets, use = "complete",
  resid.only = TRUE, ...) {

  stopifnot(class(eset) == "ExpressionSet")
  use <- match.arg(use, c("complete", "pairwise.complete"))
  
  if(missing(probesets)) {
    probesets <- featureNames(eset)
  }

  # Isolate formula terms
  lm.f <- as.formula(match.call()[["formula"]])
  lm.terms <- attr(terms(lm.f), "term.labels")
  
  # Reassemble formulate without response term
  lm.f <- reformulate(lm.terms, response = "expression")
  
  # Expression matrix
  exp.mat <- t(exprs(eset)[probesets,])
  
  if(use == "complete") {
    expression <- exp.mat
    lm.out <- lm(lm.f, data = pData(eset))
    
    if(resid.only) {
      return(t(resid(lm.out)))
    }
    
  } else if(use == "pairwise.complete") {
    lm.out <- list()
    
    for(p in probesets) {
      expression <- exp.mat[, p]
      lm.out[[p]] <- lm(lm.f, data = pData(eset), na.action = na.exclude)
    }
      
    if(resid.only) {
      r.mat <- lapply(lm.out, function(x) resid(x)[sampleNames(eset)])
      return(do.call("rbind", r.mat))
    }
  }
  
  return(lm.out)
}

