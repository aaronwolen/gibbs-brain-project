# Print html tables for wiki
#' @param display named vector indicating which display option to use for each
#'   corresponding column
print_xtable <- function(x, col.names, display) {
  require(xtable, quietly = TRUE)
  
  if(class(x) != "data.frame") {
    x <- data.frame(x)
  }
  
  x <- xtable(x)
  
  if(!missing(display)) {
    x.display <- tail(display(x), -1)
    names(x.display) <- colnames(x)
    x.display[names(display)] <- display
    display(x) <- c("s", x.display)
  }
  
  if(!missing(col.names)) {
    if(is.logical(col.names)) {
      add.colnames <- col.names
    } else {
      colnames(x) <- col.names 
    }
  } else {
    add.colnames <- FALSE
  }
  
  print(x, type = "html", 
        include.colnames = ifelse(add.colnames, TRUE, FALSE))
}
