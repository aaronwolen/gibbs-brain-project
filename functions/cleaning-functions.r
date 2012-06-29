# grab my dataframetools gist (function requires devtools package)
devtools:::source_gist("https://raw.github.com/gist/540b67c86c3c81220462/785b87f8817cff500cd9d5be1c51a082db1cc842/dataframetools.r")


#' Clean-up annoying format of phenoData characteristic columns
#' 
#' The provided phenoData contains several columns with useful phenotypic
#' information such as age and tissue. These have unhelpful column names like
#' "characteristics_ch1.1" and the values are formatted so that each cell
#' contains the phenotype identifier and associated value (e.g., gender: male).
#' This function reformats these columns by keeping only the values following
#' the colon and uses the identifier preceeding the colon as the new column
#' name. This is only done if values preceeding the colons are invariant within
#' a column.
#' 
#' @param df data.frame
reformat_charc_cols <- function(df) {
  
  # Identify characteristic columns
  char.cols <- grep("characteristic", colnames(df))
  
  # Split char column values by colons
  parsed.cols <- sapply(char.cols, simplify = F, function(x) 
    do.call("rbind", strsplit(as.character(df[, x]), "\\s*:\\s*")))
  
  # Only reformat columns if values left of colon are invariant
  keep <-unlist(lapply(parsed.cols, function(x)
    length(unique(x[, 1])) == 1))
  
  char.cols <- char.cols[keep]
  parsed.cols <- parsed.cols[keep]
  
  # Loop through each parsed char column and replace corresponding column 
  # in original df with values right of colon and rename using invariant value
  # left of colon
  for(i in 1:length(char.cols)) {
    
    orig <- names(df)[char.cols[i]] 
    
    df[, char.cols[i]] <- parsed.cols[[i]][, 2]
    names(df)[char.cols[i]] <- parsed.cols[[i]][1, 1]
    
    message(paste(orig, "renamed to", names(df)[char.cols[i]]))
    
  }
  
  return(df)
}



#' Parse genomic coordinates into separate variables
#' 
#' @param x character vector of genomic coordinates to be parsed
#' @param pattern character string containing pattern used to store coordinates
#'   as a single string
#'  
#' @example parse_coords("16:28797486-28797825")

parse_coords <- function(x, pattern = "chr:start-end") {
  
  x <- as.character(x)
  
  splits <- paste(
    unlist(strsplit(gsub("\\d|\\w", "", pattern), "")), 
    collapse = "|")
  
  labels <- unlist(strsplit(pattern, split = "\\W"))
  
  parsed <- data.frame(do.call("rbind", strsplit(x, split = splits)),
                       stringsAsFactors = FALSE, row.names = NULL)
  names(parsed) <- labels
  
  if(grepl("start", pattern)) {
    parsed$start <- as.numeric(parsed$start)
  }
  
  if(grepl("end", pattern)) {
    parsed$end <- as.numeric(parsed$end)
  }
  
  return(parsed)
}



#' Grab first annotation value from multiple delimited values
#' 
#' @param x character vector
#' @param sep delimiting character
#' @param ... additional arguments passed to strsplit
#'
 
grab_first <- function(x, sep = ":", ...) {

 multi.rows <- grep(sep, x)
 
 split.rows <- strsplit(x[multi.rows], split = sep, ...)
 split.rows <- unlist(lapply(split.rows, function(y) y[1]))
 
 x[multi.rows] <- split.rows
 return(x)
}
