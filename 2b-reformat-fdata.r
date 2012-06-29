# Clean-up feature annotation data 


# Tidy expression fdata ---------------------------------------------------

fdata.exp <- fData(geo.exp[[1]])
fmeta.exp <- fvarMetadata(geo.exp[[1]])

# Change metadata to match recommended format
fmeta.exp <- fmeta.exp[, c("Column", "Description")]
fmeta.exp <- rename(fmeta.exp, 
  c("Column" = "label", "Description" = "labelDescription"))

# Invariant columns
head(fdata.exp[,invariant_cols(fdata.exp)])
fdata.exp <- fdata.exp[, !names(fdata.exp) %in% invariant_cols(fdata.exp)]

# Redundant columns
redundant_cols(fdata.exp)
fdata.exp <- subset(fdata.exp, 
  select = -c(Source_Reference_ID, Accession, GB_ACC, ILMN_Gene))

# Update metadata
fmeta.exp <- subset(fmeta.exp, label %in% names(fdata.exp))

# Add my standard fdata columns (or rename if they already exist)
# (symbol, chr, strand, position)
exp.renames <- c("Symbol" = "symbol", "Chromosome" = "chr", 
                 "Probe_Chr_Orientation" = "strand")

fdata.exp <- rename(fdata.exp, exp.renames)
fmeta.exp$label[match(names(exp.renames), fmeta.exp$label)] <- exp.renames

# Extract position from 'Probe_Coordinates'
fdata.exp$position <- gsub("(\\d)-.*", "\\1", fdata.exp$Probe_Coordinates)
fmeta.exp <- rbind(fmeta.exp, data.frame(label = "position", 
  labelDescription = "Chromosome probe position"))

fdata.exp <- reorder_cols(fdata.exp, 
  c("symbol" = 2, "chr" = 3, "position" = 4, "strand" = 5))
fmeta.exp <- fmeta.exp[match(names(fdata.exp), fmeta.exp$label),]

# Tidy column names by converting to lower case
names(fdata.exp) <- tolower(names(fdata.exp))
fmeta.exp$label <- tolower(fmeta.exp$label)

# Properly classify columns
fdata.exp$chr <- factor(fdata.exp$chr, levels = c(1:22, "X", "Y"))
fdata.exp$position <- as.numeric(as.character(fdata.exp$position))

# Rename metadata columns
rownames(fmeta.exp) <- fmeta.exp$label



# Tiday methylation fdata -------------------------------------------------

fdata.meth <- fData(geo.meth[[1]])
fmeta.meth <- fvarMetadata(geo.meth[[1]])

# Change metadata to match recommended format
fmeta.meth <- fmeta.meth[, c("Column", "Description")]
fmeta.meth <- rename(fmeta.meth, 
  c("Column" = "label", "Description" = "labelDescription"))

# Invariant columns
head(fdata.meth[,invariant_cols(fdata.meth)])
fdata.meth <- fdata.meth[, !names(fdata.meth) %in% invariant_cols(fdata.meth)]

# Redundant columns
redundant_cols(fdata.meth)
fdata.meth <- subset(fdata.meth, select = -c(GB_ACC, Name))

# Update metadata
fmeta.meth <- subset(fmeta.meth, label %in% names(fdata.meth))

# Add my standard fdata columns (or rename if they already exist)
# (symbol, chr, strand, position)
meth.renames <- c("Symbol" = "symbol", "Chr" = "chr", 
                 "RANGE_STRAND" = "strand", "MapInfo" = "position")

fdata.meth <- rename(fdata.meth, meth.renames)
fmeta.meth$label[match(names(meth.renames), fmeta.meth$label)] <- meth.renames

fdata.meth <- reorder_cols(fdata.meth, 
  c("symbol" = 2, "chr" = 3, "position" = 4, "strand" = 5))
fmeta.meth <- fmeta.meth[match(names(fdata.meth), fmeta.meth$label),]

# Tidy column names by converting to lower case
names(fdata.meth) <- tolower(names(fdata.meth))
fmeta.meth$label <- tolower(fmeta.meth$label)

# Properly classify columns
fdata.meth$chr <- factor(fdata.meth$chr, levels = c(1:22, "X", "Y"))
fdata.meth$position <- as.numeric(as.character(fdata.meth$position))

# Rename metadata columns
rownames(fmeta.meth) <- fmeta.meth$label


# Tidy microRNA fdata ----------------------------------------------------

fdata.micro <- defactor(fData(geo.micro[[1]]))
fmeta.micro <- fvarMetadata(geo.micro[[1]])

# Change metadata to match recommended format
fmeta.micro <- fmeta.micro[, c("Column", "Description")]
fmeta.micro <- rename(fmeta.micro, 
  c("Column" = "label", "Description" = "labelDescription"))

# Invariant columns
head(fdata.micro[,invariant_cols(fdata.micro)])

# Redundant columns
redundant_cols(fdata.micro)
fdata.micro <- subset(fdata.micro, 
  select = -c(Search_Key, TargetMatureName, SYMBOL))

# Update metadata
fmeta.micro <- subset(fmeta.micro, label %in% names(fdata.micro))

# Add my standard fdata columns (or rename if they already exist)
# (symbol, chr, strand, position)
micro.renames <- c("ILMN_Gene" = "symbol")

fdata.micro <- rename(fdata.micro, micro.renames)
fmeta.micro$label[match(names(micro.renames), fmeta.micro$label)] <- micro.renames

# Extract first chr from available data
fdata.micro$chr <- grab_first(fdata.micro$Chromosome, ",")
fdata.micro$strand <- grab_first(fdata.micro$Probe_Chr_Orientation, ",")

# Extract first position from available data
    
# Ridiculously, multiple microRNA positions are comma delimited in addition to
# the numbers being comma-formatted. For example, I'm assuming that 
# "100,576,207,100,576,000" represents 100576207 and 100576000
  "100,562,923,100,563,000"
parse_coords_micro <- function(x) {
  
  hits <- grep(",", x)
  parsed <- strsplit(x[hits], split = ",")
  
  make_number <- function(y) {
    # Number isn't comma formatted
    if(any(nchar(y) > 3)) return(y[1]) 
    # Must represent a number less than 1 million since no chr is > 350 Mb
    if(as.numeric(y[1]) > 350) return(paste(y[1:2], collapse = ""))
    # Else paste first three groups of three numbers to together
    return(paste(y[1:3], collapse = ""))
  }
  
  parsed <- unlist(lapply(parsed, make_number))
  
  x[hits] <- unlist(parsed)
  return(as.numeric(x))
}  

# Parse the microRNA positions and perform sanity check
pos.micro <- parse_coords_micro(fdata.micro$Probe_Coordinates)

if(any(na.omit(pos.micro) > 3.5 * 10^8)) {
  stop("parsed_coord_micro produced microRNA positions > 350 Mb")
}

fdata.micro$position <- pos.micro

# Add new columns to metadata
fmeta.micro <- rbind(fmeta.micro, data.frame(
  label = c("chr", "strand", "position"), 
  labelDescription = paste("First available microRNA", 
    c("chromosome", "chromosome position", "strand"))))

fdata.micro <- reorder_cols(fdata.micro, 
  c("symbol" = 2, "chr" = 3, "position" = 4, "strand" = 5))
fmeta.micro <- fmeta.micro[match(names(fdata.micro), fmeta.micro$label),]

# Tidy column names by converting to lower case
names(fdata.micro) <- tolower(names(fdata.micro))
fmeta.micro$label <- tolower(fmeta.micro$label)

# Properly classify columns
fdata.micro$chr <- factor(fdata.micro$chr, levels = c(1:22, "X", "Y"))
fdata.micro$position <- as.numeric(as.character(fdata.micro$position))

# Rename metadata columns
rownames(fmeta.micro) <- fmeta.micro$label

