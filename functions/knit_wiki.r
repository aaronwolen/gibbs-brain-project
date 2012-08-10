
#' @param filename input filename
#' @param output.path output file path
#' @param output file type 
#' @param output.filename (optional) different basename for output file. This
#'   may also be a function that manipulates the input filename

knit_wiki <- function(filename, output.path = "wiki", output.ext = "md",
  output.filename) {
  require(knitr, quietly = TRUE)
  
  # Create output filename and add proper extension 
  output <- ifelse(missing(output.filename), filename, output.filename)
  output <- gsub("\\.\\w+$|$", paste("\\.", output.ext, sep = ""), output)
  
  # knit document and save output in wiki directory
  output <- file.path(output.path, output)
  
  
  knit(filename, output)
}
