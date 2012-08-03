# Extract phenodata from each of the geo datasets and combine into a single 
# data.frame containing sample phenotypic data to be added to each expressionSet
# object's pData slot.

# Use appropriate columns from pData to constrcut MIAME objects

# Construct common phenotype data.frame -----------------------------------

pdata <- do.call("rbind", lapply(geo.data, pData))
pdata <- defactor(pdata)

# Reformat annoyingly formatted 'characteristic' columns
pdata <- reformat_charc_cols(pdata)


# Create identifiers for individuals and tissue samples -------------------

# Parse out individual identifiers from 'title'
id.cols <- data.frame(do.call("rbind", strsplit(pdata$title, "-")))
names(id.cols) <- c("tissuebank", "individual", "tissue", "assay")

# merge parsed identifiers from 'title' after removing the now redundant
# 'tissue' and 'tissuebank' columns
pdata <- subset(pdata, select = -c(tissue, tissuebank))
pdata <- cbind(id.cols, pdata)

# remove the now unnecessary title column
pdata <- subset(pdata, select = -title)

# Ensure geo_accession column provides unique identifiers
nrow(pdata)  == length(unique(pdata$geo_accession))
rownames(pdata) <- pdata$geo_accession

# Ensure 'individual' column is a unique identifier within each assay/tissue
check_max.n <- nlevels(pdata$tissue) * nlevels(pdata$assay)
if(any(!table(pdata$individual) <= check_max.n)) {
  stop("pdata$individual is not a unique identifier.")
}


# Isolate invariant and MIAME-relevant columns ----------------------------
invar.pdata <- pdata[1, invariant_cols(pdata)]

# Isolate other MIAME-relevant columns
miame.cols <- c("submission_date", "molecule_ch1", "extract_protocol_ch1", 
                "label_ch1", "label_protocol_ch1", "hyb_protocol", "scan_protocol", 
                "data_processing", "platform_id")
miame.pdata <- pdata[, miame.cols]

pdata <- pdata[, setdiff(names(pdata), c(names(invar.pdata), miame.cols))]

# Remove columns redundantly identifying tissue type or source molecule
with(pdata, cbind(source_name_ch1, description, tissue))

pdata <- subset(pdata, select = -c(type, description, source_name_ch1))

# Remove unnecessary data_row_count column
pdata <- subset(pdata, select = -data_row_count)


#  Add in supplementary table 1 data --------------------------------------

sample.pdata$individual <- gsub("\\D", "", sample.pdata$sample)

# Any missing samples?
setdiff(sample.pdata$individual, pdata$individual)
setdiff(pdata$individual, sample.pdata$individual)

missing.inds <- which(!sample.pdata$individual %in% pdata$individual)

pdata <- merge(pdata, sample.pdata[, c("individual", "death")], by = "individual")


# Rename and reorder columns ----------------------------------------------
pdata <- rename(pdata, c("geo_accession" = "geo", "prep_hyb_batch" = "batch"))

order.cols <- c("geo", "individual", "assay", "tissue", "tissuebank",
  "batch", "gender", "age", "pmi", "death")

if(any(!colnames(pdata) %in% order.cols)) {
  stop("Columns will be lost after reordering.")
}

pdata <- pdata[, order.cols]


# Properly set column classes ---------------------------------------------
pdata <- classify_columns(pdata, 
  num.cols = c("age", "pmi"), 
  fac.cols = c("tissue", "gender", "individual", "tissuebank", "batch"))

pdata$gender <- factor(pdata$gender, levels = c("male", "female"))
pdata$individual <- factor(pdata$individual, 
  levels = sort(as.numeric(levels(pdata$individual))))

# Reorder
pdata <- arrange(pdata, assay, tissue, individual)

# Rownames of pdata must match colnames of expression matrix
rownames(pdata) <- pdata$geo

write.csv(pdata, "data/cleaned-pdata.csv")


# Assay-specific pdata files ----------------------------------------------

pdata.exp <- subset(pdata, assay == "mRNA")
pdata.meth <- subset(pdata, assay == "CpG")
pdata.micro <- subset(pdata, assay == "miRNA")

# Construct general MIAME object ------------------------------------------

# Use data in invar.pdata
miame <- MIAME(
  name = invar.pdata$contact_name,
   lab = invar.pdata$contact_laboratory,
   contact = invar.pdata$contact_email,
   title = "Abundant Quantitative Trait Loci for CpG Methylation and Expression Across Human Brain Tissues", 
   url = invar.pdata$contact_web_link,
   samples = as.list(pdata$geo),
   pubMedIds = "20485568",
   other = list(status = invar.pdata$status,
                updated = invar.pdata$last_update_date))

miame@name <- gsub(",+", " ", miame@name)


# Assay-specific MIAME objects --------------------------------------------

# This function will place the miame.pdata columns in their proper slots
fortify_miame <- function(miame, pdata) {
  miame@hybridizations <- list(extract = pdata$extract_protocol_ch1[1],
                         label = pdata$label_ch1[1],
                         scan = pdata$scan_protocol[1])
   miame@preprocessing <- list(protocol = pdata$label_protocol_ch1[1],
                        processing = pdata$data_processing[1])
  miame@other = c(miame@other, list(submission = pdata$submission_date[1]))
  
  return(miame)
}

miame.exp <- fortify_miame(miame, subset(miame.pdata, platform_id == gpl["exp"]))
miame.meth <- fortify_miame(miame, subset(miame.pdata, platform_id == gpl["meth"]))
miame.micro <- fortify_miame(miame, subset(miame.pdata, platform_id == gpl["micro"]))
