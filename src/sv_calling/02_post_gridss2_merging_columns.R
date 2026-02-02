library(stringr)

args <- commandArgs(trailingOnly = TRUE)
input_path <- args

# Step 1: Read the VCF file and separate header comments
print(paste0("Processing sample ", input_path))
lines <- readLines(input_path)
header_lines <- lines[grep("^##", lines)]
content_lines <- lines[!grepl("^##", lines)]

# Convert content lines into a data frame
df_content <- read.csv(textConnection(content_lines), header = FALSE, sep = "\t", as.is = TRUE, comment.char = "#")

# Step 2: Use a more efficient approach to merge SV data from tumor samples
calculate_pbk_value <- function(infoline_value, df_row, numcol) {
  column_indexes <- 11:numcol
  
  if (infoline_value %in% c('SR', 'RP', 'REF', 'REFPAIR', 'BASRP', 'BASSR', 'BANRP', 'BANSR', 'BSC', 'BUM', 'VF', 'BVF', 'ASRP', 'ASSR')) {
    value2 <- sum(sapply(column_indexes, function(k) as.numeric(strsplit(df_row[k], ":", fixed = TRUE)[[1]][which(infoline == infoline_value)])))
    return(value2)
  } else {
    value2 <- max(sapply(column_indexes, function(k) as.numeric(strsplit(df_row[k], ":", fixed = TRUE)[[1]][which(infoline == infoline_value)])), na.rm = TRUE)
    return(value2)
  }
}

infoline <- strsplit(df_content$V9[1], ":", fixed = TRUE)[[1]]
numcol <- ncol(df_content)

df_content$pbk <- apply(df_content, 1, function(df_row) {
  pbk_values <- sapply(infoline[-1], function(infoline_value) calculate_pbk_value(infoline_value, df_row, numcol))
  return(paste0(".", ":", paste0(pbk_values, collapse = ":")))
})

df_content <- df_content[, c(paste0("V", c(1:10)), "pbk")]
colnames(df_content) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "COMPOSITE_TUMOR_LIBRARIES")

# Step 3: Write the output file, preserving the original header and appending the processed content
output_path <- gsub("\\.vcf.gz$", ".composite.v5.vcf.gz", input_path)
gz_output <- gzfile(output_path, "w")
writeLines(header_lines, gz_output)
write.table(df_content, file = gz_output, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE)
close(gz_output)

