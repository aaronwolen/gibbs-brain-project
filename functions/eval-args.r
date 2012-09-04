# Evaluate supplied command line arguments
args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0){
  stop("No supplied arguments.")
} else {
  for(i in seq_along(args)) {
    eval(parse(text = args[[i]]))
  }
}
  