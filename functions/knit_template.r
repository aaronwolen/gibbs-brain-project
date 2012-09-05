#' Use brew to programmatically set knitr chunk labels 

#' @param file a character string naming the file to be read in. 
#' @param label character string matching text within brew delimiters
#' @param labels character vector containing values to assign to label

knit_template <- function(file, variable = "variable", labels) {
  require(knitr, quietly = TRUE)

  # Set knitr to read brew syntax
  pat_brew()
  knit_hooks$set(inline = identity)
    
  # Process template file for each level of supplied labels
  src <- NULL
  for(i in labels) {
    assign(variable, i)
    src <- c(src, knit(text = readLines(file)))
  }
  
  # Switch back to markdown syntax
  pat_md()
  render_markdown()
  out <- knit(text = src)
  
  return(out)
}
