# Print html tables for wiki
#' @param display named vector indicating which display option to use for each
#'   corresponding column
#'   
print_xtable <- function(x, display, col.names = TRUE, row.names = FALSE) {
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
  
  if(!missing(row.names)) {
    if(is.logical(row.names)) {
      add.rownames <- row.names
    } else {
      rownames(x) <- row.names
    }
  } else {
    add.rownames <- FALSE
  }
  
  out <- print.xtable(x, type = "html", 
        include.colnames = add.colnames,
        include.rownames = add.rownames,
        print.results = FALSE)
  
  # Strip out html comments and trailing newline
  out <- gsub("<!.*-->", "", out)
  out <- gsub("^\n|\n$", "", out)
  
  # Add missing semicolon to html symbols
  out <- gsub("(&\\w+)", "\\1;", out)
  
  cat(out)
}
