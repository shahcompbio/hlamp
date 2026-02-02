library(dplyr)
library(stringr)
library(signals)
library(data.table)
library(anndata)
library(scales)
library(viridisLite)

args <- commandArgs(trailingOnly = TRUE)
w <- args

print(paste0("Processing ", w))

options(scipen = 999)
hg19 <- read.csv("./references/hg19.coordinates.csv")
blklst <- read.csv("./references/encode_dac_blacklist/ENCFF001TDO.nochr.bed", header=F, as.is=T, sep="\t")
blklst <- blklst[,c(1:4)]
colnames(blklst) <- c("chr", "start", "end", "class")
qc <- fread(paste0("./analysis/metadata/qc/", w, "/", w, "_qc.csv.gz"))
qc$cell_id <- str_sub(qc$cell_id, -15)
ad <- read_h5ad(paste0("./sv_gc_normalized/", w, "/data/", w, "_corrected.h5ad"))
cl <- fread(paste0("./sv_gc_normalized/", w, "/tables/", w, "_clusters_0.7_1000.csv"))
cl$cell_id <- str_sub(cl$cell_id, -15)
if (length(cl$cell_id) != length(unique(cl$cell_id))){
  stop("Alert: redundant cell barcodes.")
}

modcol <- c("#3182BD","#9ECAE1","#CCCCCC","#FDCC8A","#F79763","#F2663D","#E42127","#BB1E3A","#961C46","#6A2539","#412028","#000000")
setkey(cl, cluster_label)

# Assuming hg19 is a data.frame, convert it to a data.table
hg19_dt <- as.data.table(hg19)

# Define the chromosome sets
autosomes_long_plus_X <- as.character(c(1:12, 16:20, "X"))
autosomes_acro <- as.character(c(13:15, 21, 22))

# Extract the relevant centromere coordinates for merging
# We only need the start/end coordinates once per chromosome
hg19_centro_coords <- hg19_dt[, .(chr, centro_start, centro_end)]
hg19_centro_coords$centro_start <- hg19_centro_coords$centro_start - 1000000
hg19_centro_coords$centro_end <- hg19_centro_coords$centro_end + 1000000

# Chromosomes with p arm artifact
hg19_centro_coords$centro_start[hg19_centro_coords$chr == "2"] <- hg19_centro_coords$centro_start[hg19_centro_coords$chr == "2"] - 3000000
hg19_centro_coords$centro_start[hg19_centro_coords$chr == "16"] <- hg19_centro_coords$centro_start[hg19_centro_coords$chr == "16"] - 3000000
hg19_centro_coords$centro_start[hg19_centro_coords$chr == "9"] <- hg19_centro_coords$centro_start[hg19_centro_coords$chr == "9"] - 8000000

# Chromosomes with q arm artifact
hg19_centro_coords$centro_end[hg19_centro_coords$chr == "1"] <- hg19_centro_coords$centro_end[hg19_centro_coords$chr == "1"] + 25000000
hg19_centro_coords$centro_end[hg19_centro_coords$chr == "9"] <- hg19_centro_coords$centro_end[hg19_centro_coords$chr == "9"] + 18000000
hg19_centro_coords$centro_end[hg19_centro_coords$chr == "16"] <- hg19_centro_coords$centro_end[hg19_centro_coords$chr == "16"] + 8000000

# Extract encode blacklist regions
blklst <- as.data.table(blklst)
encode_blacklist_coords <- blklst[, .(chr, start, end)]

# Extra filter
xfilter_coords <- as.data.frame(matrix(c("4", "17", 68260000, 41380000, 68270000, 41480000), nrow = 2, ncol = 3)) # these two intervals are likely associated with unequal coverage of the unmatched normal cells.
colnames(xfilter_coords) <- c("chr", "start", "end")
xfilter_coords$start <- as.numeric(xfilter_coords$start)
xfilter_coords$end <- as.numeric(xfilter_coords$end)
xfilter_coords <- as.data.table(xfilter_coords)

# Define the set of 'hq tumor cells' (or whichever subset you need)
target_cell_ids <- qc %>%
  filter(
    quality >= 0.75,
    total_mapped_reads > 250000,
    is_control == FALSE,
    is_contaminated == FALSE,
    is_tumor_cell == "yes",
    is_hq == "yes"
  ) %>%
  pull(cell_id)

# -------------------------------------------------------------
# STEP 1: Bulk Extraction of Data
# -------------------------------------------------------------

# Subsetting the AnnData object to ONLY include the target cells (ROWS)
# This is much faster than iterating and subsetting inside a loop.
ad_subset <- ad[target_cell_ids[target_cell_ids %in% row.names(ad$obs)], ]

# Extract ALL layers/data for the target cells in one go
# AnnData objects are (Cell x Bins), we need to TRANSPOSE them (t())
# so that the resulting matrix columns are the BINS, and we can bind them
# to the ad$var (Bin-wise metadata) later.

# The extracted layer data is often a sparse matrix, which is efficient.
# Use as.matrix() to convert them to standard dense matrices for data.table/dplyr.
cn_matrix_t <- t(as.matrix(ad_subset$layers["gc_corrected"]))
state_matrix_t <- t(as.matrix(ad_subset$layers["rounded_gc_corrected"]))
normal_reads_matrix_t <- t(as.matrix(ad_subset$layers["reads_normal"]))
tumor_reads_matrix_t <- t(as.matrix(ad_subset$X))

# -------------------------------------------------------------
# STEP 2: Convert to a Long Format Data Frame (The Main Efficiency Gain)
# -------------------------------------------------------------

# Get the genomic bin information (the "var" part)
df_var <- as.data.table(ad_subset$var)
df_var[, bin_id := rownames(ad_subset$var)] # Add the bin name/index

# Get the cell ID information (the "obs" part)
cell_ids <- rownames(ad_subset)
id_cols <- names(df_var)

# Use data.table to efficiently combine and melt the data
# 1. Create a base data.table from the genomic metadata (df_var)
# 2. Add the CN matrix, which is aligned with the rows of df_var
# 3. Melt the data to convert from a wide (bins as columns) to a long (one row per bin-cell pair) format

df <- data.table(
  df_var,
  copy = cn_matrix_t,
  state = state_matrix_t,
  normal_reads = normal_reads_matrix_t,
  reads = tumor_reads_matrix_t
)

df_long <- melt(df,
                    # ID columns to keep (genomic metadata)
                    id.vars = id_cols,
                    
                    # Columns to melt (all the copy/state/reads columns)
                    # If you have non-ID, non-CN columns, you would list them explicitly.
                    # Since you only added CN columns, we can skip 'measure.vars'
                    
                    # Tell melt how to separate the old column name into 'variable' (cell_id) and 'layer' (copy/state/reads)
                    # The pattern says: split the column name into two parts (e.g., 'copy' and '130057A-...')
                    measure.vars = patterns("copy\\.|state\\.|reads\\.|normal_reads\\."),
                    
                    variable.name = "cell_id_index",  # Temporary: will hold the index of the cell
                    value.name = "value_placeholder"  # Temporary: will hold the data value
)

setnames(df_long, "cell_id_index", "full_col_name")

# 1. SPLIT: Separate the Layer Name (e.g., 'copy') from the Cell ID (e.g., '130057A...')
# We use the 'tstrsplit' function which is data.table's fast string splitter.
df_long[, c("layer_name", "cell_id") := tstrsplit(full_col_name, ".", fixed = TRUE)]

# 2. DROP: Remove the columns we no longer need.
df_long[, c("full_col_name") := NULL]

# 3. PIVOT: Use dcast to pivot the 'layer_name' column to become the new value columns.
# The row-defining columns are ALL the genomic metadata (id_cols) PLUS the new 'cell_id'.
# The column-defining variable is 'layer_name' (copy, state, normal_reads).
df_final <- dcast(df_long,
                      formula = ... + cell_id ~ layer_name,
                      value.var = "value_placeholder"
)


# 4. Filter and select the final columns to match your desired output

# Assuming 'hg19' is available, run your vectorized operations
df_final[, `:=`(
    plot.start = start + hg19$add_values[match(chr, hg19$chr)],
    plot.end = end + hg19$add_values[match(chr, hg19$chr)],
    segment_length = end - start + 1
)]
# (You can re-add this block when you have the hg19 object available)


# 5. Vectorized Update for centromeres and acrocentric arms (All Logic in One Block)
df_final[, centroacro := "no"]
df_final[
  hg19_centro_coords,
  on = "chr",
  centroacro := fifelse(
    (chr %in% autosomes_long_plus_X & x.start <= i.centro_end + 1000000 & x.end >= i.centro_start - 1000000) |
      (chr %in% autosomes_acro & x.start <= i.centro_end),
    "yes",
    # x.centroacro refers to the existing value of centroacro in df_final
    x.centroacro
  )
]

df_final[, blacklist := "no"]
df_final[
  encode_blacklist_coords,
  on = "chr",
  blacklist := fifelse(
    (chr %in% c(as.character(c(1:22)), "X", "Y") & x.start <= i.end & x.end >= i.start),
    "yes",
    x.blacklist
  )
]

df_final[, xfilter := "no"]
df_final[
  xfilter_coords,
  on = "chr",
  blacklist := fifelse(
    (chr %in% c("4", "17") & x.start <= i.end & x.end >= i.start),
    "yes",
    x.xfilter
  )
]

# Select and order the final columns
df_final <- df_final[, .(
  chr,
  start,
  end,
  reads,
  copy,
  state,
  cell_id,
  normal_reads,
  bin_id,
  segment_length,
  plot.start,
  plot.end,
  centroacro,
  blacklist,
  xfilter
)]
df_final[, sort_rank := match(cell_id, cl$cell_id)]
setorder(df_final, sort_rank, chr, start)

fwrite(df_final, paste0("./analysis/normalization/cn_dataframe/", w, ".normalized.gc_corrected.filters.annotated.csv"), quote = F, row.names = F, col.names = T, sep = ",")

# Apply custom filtering
df_final <- df_final[
  normal_reads != 0 &
    segment_length >= 1000 &
    centroacro == "no" &
    blacklist == "no" &
    xfilter == "no"
]

# Annotating clone assignment data
infoline <- c("cell_id", "sample_id", "sample_type", "quality", "coverage_depth", "coverage_breadth")
qc_subset <- qc[, ..infoline]
cl[qc_subset,
   # := creates the new columns in the cl data table
   `:=`(
     sample_id = i.sample_id,
     sample_type = i.sample_type,
     quality = i.quality,
     coverage_depth = i.coverage_depth,
     coverage_breadth = i.coverage_breadth
   ),
   # on specifies the column to match (no need to explicitly setkey(cl))
   on = "cell_id"
]

cl[df_final,
   `:=`(
     sort_rank = i.sort_rank
   ),
   # on specifies the column to match (no need to explicitly setkey(cl))
   on = "cell_id"
]

cl[, quality_color := col_numeric(
  palette = viridisLite::cividis(256), 
  domain = c(0.75, 1.0) # <--- This one change scales the colors as requested
)(quality)]

# Consensus copy number
tempdf <- consensuscopynumber(as.data.frame(df_final[,c("chr", "start", "end", "reads", "copy", "state", "cell_id")]))
tempdf$plot.start <- 0
tempdf$plot.end <- 0
for (i in unique(tempdf$chr)){
  tempdf$plot.start[tempdf$chr == i] <- tempdf$start[tempdf$chr == i] + hg19$add_values[hg19$chr == i]
  tempdf$plot.end[tempdf$chr == i] <- tempdf$end[tempdf$chr == i] + hg19$add_values[hg19$chr == i]
}

for (i in unique(cl$cluster_label)){
  cscn <- consensuscopynumber(as.data.frame(df_final[cell_id %in% cl$cell_id[cl$cluster_label == i],c("chr", "start", "end", "reads", "copy", "state", "cell_id")]))
  tempdf[,paste0("clone_", i)] <- cscn$copy
  rm(cscn)
}

tempdf$length <- tempdf$end - tempdf$start
tempdf <- tempdf[tempdf$length >= 1000,]
write.csv(tempdf, paste0("./analysis/normalization/subclonal_cn_v4/", w, "_clusters_0.7_1000_subclonal_cn.gc_corrected.csv"), quote = F, row.names = F)
#old_par <- par(no.readonly = TRUE)

pdf(paste0("./figures/composite_plot_dlp_v4/composite.dlp.pseudobulk.sc.sv-seg-cn.gc_corrected.", w, ".pdf"), height=10, width=12)
layout_matrix <- matrix(c(1, 2, 3), ncol = 1) 
# Heights: 35% for Plot 1, 65% for Plot 2
layout(layout_matrix, heights = c(0.25, 0.5, 0.25))

par(mar = c(0.5, 4.1, 4.1, 2.1))
plot(NULL, xlim = c(hg19$add_values[1], hg19$add_values[nrow(hg19)]+hg19$length[nrow(hg19)]), ylim = c(0, 30), ylab = "Copy number", xlab = "", las=1, xaxt = "n", frame = FALSE, main = paste0("Normalized, GC-corrected, Filtered CN in ", w, "; # hq cancer cells = ", max(df_final$sort_rank), "\nheatmap -- clone;sample_id;sample_type;quality(0.75-1.00)"))
axis(side = 1, labels = F, at = hg19$add_values)
#text(x = hg19$text_location[1:24], y = -3, labels = hg19$chr[1:24], xpd=NA)
abline(v = hg19$add_values, col = "gray")
tempdf$state[tempdf$state < 0] <- 0
tempdf$state[tempdf$state >= 11] <- 11
segments(x0 = tempdf$plot.start[tempdf$copy != "NaN" & tempdf$copy < 30], y0 = tempdf$copy[tempdf$copy != "NaN" & tempdf$copy < 30], x1 = tempdf$plot.end[tempdf$copy != "NaN" & tempdf$copy < 30], y1 = tempdf$copy[tempdf$copy != "NaN" & tempdf$copy < 30], col = modcol[tempdf$state[tempdf$copy != "NaN" & tempdf$copy < 30]+1], lwd = 2)
if (nrow(tempdf[tempdf$copy != "NaN" & tempdf$copy >= 30,]) != 0){
  segments(x0 = tempdf$plot.start[tempdf$copy != "NaN" & tempdf$copy >= 30], y0 = 30, x1 = tempdf$plot.end[tempdf$copy != "NaN" & tempdf$copy >= 30], y1 = 30, col = "navy", lwd = 3)
}
legend("topright", legend = names(table(tempdf$state)), fill = modcol, border = NA, bty = 'n')

# Main CN plot (lower)
par(mar = c(0.5, 4.1, 0.1, 2.1)) 
cell_number <- max(df_final$sort_rank)
plot(NULL, xlim = c(hg19$add_values[1], hg19$add_values[nrow(hg19)]+hg19$length[nrow(hg19)]), ylim = c(0, cell_number), ylab = "Cancer Cell Copy Number", xlab = "", las=1, xaxt = "n", frame = FALSE)
axis(side = 1, labels = F, at = hg19$add_values)
text(x = hg19$text_location[1:24], y = -max(df_final$sort_rank)*0.07, labels = hg19$chr[1:24], xpd=NA)
abline(v = hg19$add_values, col = "gray")
df_final$state[df_final$state < 0] <- 0
df_final$state[df_final$state >= 11] <- 11
rect(df_final$plot.start, cell_number-df_final$sort_rank, df_final$plot.end, cell_number-df_final$sort_rank+1, border = NA, col = modcol[df_final$state+1])

for (i in 1:nrow(cl)){
  rect(3100000000, cell_number-cl$sort_rank[i], 3130000000, cell_number-cl$sort_rank[i]+1, border = NA, col = viridisLite::turbo(length(names(table(cl$cluster_label))))[which(names(table(cl$cluster_label)) == cl$cluster_label[i])], xpd = T) # Clone assignment
  rect(3140000000, cell_number-cl$sort_rank[i], 3170000000, cell_number-cl$sort_rank[i]+1, border = NA, col = viridisLite::cividis(length(names(table(cl$sample_id))))[which(names(table(cl$sample_id)) == cl$sample_id[i])], xpd = T) # Sample ID
  rect(3180000000, cell_number-cl$sort_rank[i], 3210000000, cell_number-cl$sort_rank[i]+1, border = NA, col = viridisLite::cividis(length(names(table(cl$sample_type))))[which(names(table(cl$sample_type)) == cl$sample_type[i])], xpd = T) # Sample type
  rect(3220000000, cell_number-cl$sort_rank[i], 3250000000, cell_number-cl$sort_rank[i]+1, border = NA, col = cl$quality_color[i], xpd = T) # Quality
}

par(mar = c(4.1, 4.1, 0.1, 2.1))
plot(NULL, xlim = c(hg19$add_values[1], hg19$add_values[nrow(hg19)]+hg19$length[nrow(hg19)]), ylim = c(0, 20), ylab = "Copy number", xlab = "", las=1, xaxt = "n", frame = FALSE)
axis(side = 1, labels = F, at = hg19$add_values)
text(x = hg19$text_location[1:24], y = -2, labels = hg19$chr[1:24], xpd=NA)
abline(v = hg19$add_values, col = "gray")

for (i in rev(unique(cl$cluster_label))){
  segments(x0 = tempdf$plot.start[tempdf[,paste0("clone_", i)] < 20], y0 = tempdf[tempdf[,paste0("clone_", i)] < 20,paste0("clone_", i)], x1 = tempdf$plot.end[tempdf[,paste0("clone_", i)] < 20], y1 = tempdf[tempdf[,paste0("clone_", i)] < 20,paste0("clone_", i)], col = viridisLite::turbo(length(names(table(cl$cluster_label))))[which(names(table(cl$cluster_label)) == i)])
  if (nrow(tempdf[tempdf[,paste0("clone_", i)] >= 20,]) != 0){
    segments(x0 = tempdf$plot.start[tempdf[,paste0("clone_", i)] >= 20], y0 = 20, x1 = tempdf$plot.end[tempdf[,paste0("clone_", i)] >= 20], y1 = 20, col = viridisLite::turbo(length(names(table(cl$cluster_label))))[which(names(table(cl$cluster_label)) == i)], lwd = 3)
  }
}

legend("topright", legend = names(table(cl$cluster_label)), fill = viridisLite::turbo(length(names(table(cl$cluster_label)))), border = NA, bty = 'n')
dev.off()
