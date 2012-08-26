library(GEOquery)

# Create data directory ---------------------------------------------------
data.dir <- "data/geo"
dir.create(data.dir)

# Download data -----------------------------------------------------------
geo.id <- "GSE15745"

geo.data <- getGEO(geo.id, destdir = data.dir)

# Export combined geo files as single rda file
save(list = "geo.data", file = file.path("data", "geo-raw-data.rda"))


# Remove individual geo files ---------------------------------------------
geo.files <- list.files(data.dir, full.names = TRUE)
file.remove(geo.files)

# Remove geo directory
file.remove(data.dir)
