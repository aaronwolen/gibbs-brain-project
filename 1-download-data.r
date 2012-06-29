library(GEOquery)

# Download data -----------------------------------------------------------

geo.id <- "GSE15745"

geo.data <- getGEO(geo.id, destdir = "data/geo")
save(list = "geo.data", file = paste(geo.id, "raw-data.rda", sep = "_"), 
     compress = TRUE)

gsm <- getGEO(system.file("data/geo/GSE15745-GPL6104_series_matrix-1.txt.gz", package = "GEOquery"))