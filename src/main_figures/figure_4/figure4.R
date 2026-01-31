# Figure 4
library(stringr)
library(circlize)
library(ComplexHeatmap)
indi <- read.csv("./metadata/dlp.cohort.summary.individual_111925.tsv", header=T, as.is=T, sep="\t")
east <- read.csv("./modeling/amplified_regions_summary_table_011026.csv")

# Panel A, top
w <- "NCI-H69"
geneofi <- "MYCN"
svsketch.dlp.pseudobulk(w, "2", geneofi)

# Panel A, bottom
celldf <- read.csv("./for-jake/NCI-H69_ordered_loadings.csv")
qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)
qc <- qc[cell_id %in% celldf$cell_id]

colnames(celldf) <- c("cell_id", "e2", "e3", "e0", "e1")
for (i in rev(c("e3", "e2", "e1", "e0"))){
  celldf <- celldf[order(celldf[,i], decreasing = T),]
}
celldf$sort_rank <- c(nrow(celldf):1)
celldf <- celldf[order(celldf$sort_rank, decreasing = F),]
setDT(celldf)

ad <- read_h5ad(paste0("./results/sv_cna_clones_nov_27NCI-H69/data/", w, "_corrected.h5ad"))
target_cell_ids <- celldf$cell_id
ad_subset <- ad[target_cell_ids[target_cell_ids %in% row.names(ad$obs)], ]

cn_matrix_t <- t(as.matrix(ad_subset$layers["gc_corrected"]))
state_matrix_t <- t(as.matrix(ad_subset$layers["rounded_gc_corrected"]))
normal_reads_matrix_t <- t(as.matrix(ad_subset$layers["reads_normal"]))
tumor_reads_matrix_t <- t(as.matrix(ad_subset$X))

df_var <- as.data.table(ad_subset$var)
df_var[, bin_id := rownames(ad_subset$var)]

cell_ids <- rownames(ad_subset)
id_cols <- names(df_var)

df <- data.table(
  df_var,
  copy = cn_matrix_t,
  state = state_matrix_t,
  normal_reads = normal_reads_matrix_t,
  reads = tumor_reads_matrix_t
)

df_long <- melt(df,
                id.vars = id_cols,
                measure.vars = patterns("copy\\.|state\\.|reads\\.|normal_reads\\."),
                variable.name = "cell_id_index",  # Temporary: will hold the index of the cell
                value.name = "value_placeholder"  # Temporary: will hold the data value
)
setnames(df_long, "cell_id_index", "full_col_name")

# 1. SPLIT: Separate the Layer Name (e.g., 'copy') from the Cell ID (e.g., '130057A...')
df_long[, c("layer_name", "cell_id") := tstrsplit(full_col_name, ".", fixed = TRUE)]

# 2. DROP: Remove the columns we no longer need.
df_long[, c("full_col_name") := NULL]

# 3. PIVOT: Use dcast to pivot the 'layer_name' column to become the new value columns.
df_final <- dcast(df_long,
                  formula = ... + cell_id ~ layer_name,
                  value.var = "value_placeholder"
)

# 4. Filter and select the final columns to match your desired output
df_final[, `:=`(
  plot.start = start + hg19$add_values[match(chr, hg19$chr)],
  plot.end = end + hg19$add_values[match(chr, hg19$chr)],
  segment_length = end - start + 1
)]

# 5. Vectorized Update for centromeres and acrocentric arms (All Logic in One Block)
df_final[, centroacro := "no"]
df_final[
  hg19_centro_coords,
  on = "chr",
  centroacro := fifelse(
    (chr %in% autosomes_long_plus_X & x.start <= i.centro_end + 1000000 & x.end >= i.centro_start - 1000000) |
      (chr %in% autosomes_acro & x.start <= i.centro_end),
    "yes",
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
rm(cn_matrix_t, state_matrix_t, normal_reads_matrix_t, tumor_reads_matrix_t, df_var, cell_ids, id_cols, df, df_long)

df_final <- df_final[
  normal_reads != 0 &
    segment_length >= 1000 &
    centroacro == "no" &
    blacklist == "no" &
    xfilter == "no"
]
df_final[celldf, sort_rank := i.sort_rank, on = "cell_id"]

df <- east[east$individual == w,]
setDT(df)
colnames(df)[2] <- "chr"
df_final <- df_final[
  df,
  on = .(chr = chr,
         start >= start, 
         end <= end),
  j = .(
    chr = x.chr,
    start = x.start,
    end = x.end,
    reads = x.reads,
    copy = x.copy,
    state = x.state,
    cell_id = x.cell_id,
    normal_reads = x.normal_reads,
    bin_id = x.bin_id,
    segment_length = x.segment_length,
    plot.start = x.plot.start,
    plot.end = x.plot.end,
    centroacro = x.centroacro,
    blacklist = x.blacklist,
    xfilter = x.xfilter,
    sort_rank = x.sort_rank,
    regs = i.id
  ),
  
  nomatch = 0,
  allow.cartesian = TRUE
]

df_final[, chr_ord := {
  x <- sub("^chr", "", chr, ignore.case = TRUE)
  v <- suppressWarnings(as.integer(x))
  v[is.na(v) & toupper(x) == "X"] <- 23
  v[is.na(v) & toupper(x) == "Y"] <- 24
  v[is.na(v) & toupper(x) %in% c("M", "MT")] <- 25 # optional
  v[is.na(v)] <- 9999 # push any non-standard names to the end
  v
}]
setorder(df_final, sort_rank, chr_ord, start, end)
df_final[, chr_ord := NULL]

df[, chr_ord := {
  x <- sub("^chr", "", chr, ignore.case = TRUE)
  v <- suppressWarnings(as.integer(x))
  v[is.na(v) & toupper(x) == "X"] <- 23
  v[is.na(v) & toupper(x) == "Y"] <- 24
  v[is.na(v) & toupper(x) %in% c("M", "MT")] <- 25 # optional
  v[is.na(v)] <- 9999 # push any non-standard names to the end
  v
}]
setorder(df, chr_ord, start, end)
df[, chr_ord := NULL]
df$seglen <- df$end - df$start + 1
gaplen <- round(sum(df$seglen) * 0.05)
df$add_value <- 0
if (nrow(df) >= 2){
  for (i in 2:nrow(df)){
    df$add_value[i] <- df$add_value[i-1] + df$seglen[i-1] + gaplen
  }
}

regs <- df$id
for (i in 1:length(regs)){
  df_final$plot.start[df_final$regs == regs[i]] <- df_final$start[df_final$regs == regs[i]] - df$start[df$id == regs[i]] + df$add_value[df$id == regs[i]]
  df_final$plot.end[df_final$regs == regs[i]] <- df_final$end[df_final$regs == regs[i]] - df$start[df$id == regs[i]] + df$add_value[df$id == regs[i]]
}

pdf("Fig4a.bottom.pdf", width = 4, height = 4)
cell_number <- max(df_final$sort_rank)
plot(NULL, xlim = c(0, df$add_value[nrow(df)] + df$seglen[nrow(df)]), ylim = c(0, cell_number), ylab = "Cancer Cell Copy Number", xlab = "", las=1, xaxt = "n", frame = FALSE)
text(x = hg19$text_location[1:24], y = -max(df_final$sort_rank)*0.07, labels = hg19$chr[1:24], xpd=NA)
df_final$state[df_final$state < 0] <- 0
df_final$logcopy <- log10(df_final$state + 1)
df_final$linearcopy <- df_final$copy
df_final$linearcopy[df_final$linearcopy > 200] <- 200

rect(
  xleft = df_final$plot.start,
  ybottom = cell_number - df_final$sort_rank,
  xright = df_final$plot.end,
  ytop = cell_number - df_final$sort_rank + 1,
  border = NA,
  # The entire color generation logic is inside 'col'
  col = viridis::viridis_pal(option = "C", direction = 1)(100)[
    cut(
      x = df_final$linearcopy,
      # Generate 101 break points for 100 bins
      breaks = seq(0, max(df_final$linearcopy, na.rm = TRUE), length.out = 101),
      labels = FALSE,
      include.lowest = TRUE
    )
  ]
)

for (i in geneofi){
  segments(x0 = genedf$start[gsub("-", ".", genedf$gene) == i] - df$start[grepl(i, df$oncogene)] + df$add_value[grepl(i, df$oncogene)], y0 = -0.05*k, x1 = genedf$end[gsub("-", ".", genedf$gene) == i] - df$start[grepl(i, df$oncogene)] + df$add_value[grepl(i, df$oncogene)], y1 = -0.05*k, lwd = 1, xpd = T, col = "darkblue")
  text(x = genedf$start[gsub("-", ".", genedf$gene) == i] - df$start[grepl(i, df$oncogene)] + df$add_value[grepl(i, df$oncogene)], y = -0.2*k, labels = i, xpd = T)
}

dev.off()

pdf("Fig4a.bottom.label.pdf", width = 4, height = 2)
image(
  z = t(matrix(1:100, nrow = 1)),
  col = viridis::viridis_pal(option = "C", direction = 1)(100),
  axes = FALSE,
  xlab = "", ylab = ""
)

axis(
  side = 1,
  at = c(0, 100, 200, 300, 400, 500)[c(0, 100, 200, 300, 400, 500) <= max(df_final$linearcopy, na.rm = TRUE)] / 
    max(df_final$linearcopy, na.rm = TRUE),
  labels = c(0, 100, 200, 300, 400, 500)[c(0, 100, 200, 300, 400, 500) <= max(df_final$linearcopy, na.rm = TRUE)],
  
  tick = TRUE,
  line = 0
)
mtext("Copy number", side = 1, line = 2)
dev.off()

# Panel A, middle
mat <- as.matrix(celldf[, c(2:5)])
col_names <- colnames(mat) 
min_val <- min(mat)
max_val <- 100

col_list <- list(
  "e0" = circlize::colorRamp2(c(min_val, max_val), c("white", "#2fa148")), # Green
  "e1" = circlize::colorRamp2(c(min_val, max_val), c("white", "#d52c28")), # Red
  "e2" = circlize::colorRamp2(c(min_val, max_val), c("white", "#2378b5")), # Blue
  "e3" = circlize::colorRamp2(c(min_val, max_val), c("white", "#f57f20"))  # Orange
)

h_list <- NULL

for (i in seq_along(col_names)) {
  nm <- col_names[i]
  
  h <- Heatmap(
    mat[, i, drop = FALSE],
    name = nm,
    col = col_list[[nm]],
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_order = 1:nrow(mat)
  )
  
  h_list <- h_list + h
}

pdf("Fig4a.middle.loading.pdf", width = 4, height = 4)
draw(h_list)
dev.off()

# Panel A, right
w <- "NCI-H69"
i <- "130050A-R54-C64" # cell 1
i <- "130050A-R51-C43" # cell 2
i <- "130050A-R42-C46" # cell 3
i <- "130050A-R31-C57" # cell 4

z <- 4
regions <- "2-15500000-16500000;2-53000000-54000000;2-77000000-78000000;2-86000000-87000000;2-135000000-141000000;2-145000000-146000000"
geneofi <- "MYCN"
svsketch.dlp.sv.seg.cn.onecell.cssv(w, i, z, regions, geneofi)

# Panel B: conceptual illustration

# Panel C, top
w <- "GBM0510"

celldf <- read.csv("./ordered_loadings.csv")
qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)
qc <- qc[cell_id %in% celldf$cell_id]

ad <- read_h5ad(paste0("./data/", w, "_corrected.h5ad"))
target_cell_ids <- celldf$cell_id
ad_subset <- ad[target_cell_ids[target_cell_ids %in% row.names(ad$obs)], ]

cn_matrix_t <- t(as.matrix(ad_subset$layers["gc_corrected"]))
state_matrix_t <- t(as.matrix(ad_subset$layers["rounded_gc_corrected"]))
normal_reads_matrix_t <- t(as.matrix(ad_subset$layers["reads_normal"]))
tumor_reads_matrix_t <- t(as.matrix(ad_subset$X))

df_var <- as.data.table(ad_subset$var)
df_var[, bin_id := rownames(ad_subset$var)]

cell_ids <- rownames(ad_subset)
id_cols <- names(df_var)

df <- data.table(
  df_var,
  copy = cn_matrix_t,
  state = state_matrix_t,
  normal_reads = normal_reads_matrix_t,
  reads = tumor_reads_matrix_t
)

df_long <- melt(df,
                id.vars = id_cols,
                measure.vars = patterns("copy\\.|state\\.|reads\\.|normal_reads\\."),
                variable.name = "cell_id_index",  # Temporary: will hold the index of the cell
                value.name = "value_placeholder"  # Temporary: will hold the data value
)
setnames(df_long, "cell_id_index", "full_col_name")

df_long[, c("layer_name", "cell_id") := tstrsplit(full_col_name, ".", fixed = TRUE)]

df_long[, c("full_col_name") := NULL]

df_final <- dcast(df_long,
                  formula = ... + cell_id ~ layer_name,
                  value.var = "value_placeholder"
)

df_final[, `:=`(
  plot.start = start + hg19$add_values[match(chr, hg19$chr)],
  plot.end = end + hg19$add_values[match(chr, hg19$chr)],
  segment_length = end - start + 1
)]

df_final[, centroacro := "no"]
df_final[
  hg19_centro_coords,
  on = "chr",
  centroacro := fifelse(
    (chr %in% autosomes_long_plus_X & x.start <= i.centro_end + 1000000 & x.end >= i.centro_start - 1000000) |
      (chr %in% autosomes_acro & x.start <= i.centro_end),
    "yes",
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
rm(cn_matrix_t, state_matrix_t, normal_reads_matrix_t, tumor_reads_matrix_t, df_var, cell_ids, id_cols, df, df_long)

df_final <- df_final[
  normal_reads != 0 &
    segment_length >= 1000 &
    centroacro == "no" &
    blacklist == "no" &
    xfilter == "no"
]

tempdf <- consensuscopynumber(as.data.frame(df_final[,c("chr", "start", "end", "reads", "copy", "state", "cell_id")]))
tempdf$plot.start <- 0
tempdf$plot.end <- 0
for (i in unique(tempdf$chr)){
  tempdf$plot.start[tempdf$chr == i] <- tempdf$start[tempdf$chr == i] + hg19$add_values[hg19$chr == i]
  tempdf$plot.end[tempdf$chr == i] <- tempdf$end[tempdf$chr == i] + hg19$add_values[hg19$chr == i]
}
tempdf$length <- tempdf$end - tempdf$start
tempdf <- tempdf[tempdf$length >= 1000,]

pdf(paste0("Fig4c.top.", w, ".pdf"), width = 4, height = 3)
plot(NULL, xlim = c(hg19$add_values[1], hg19$add_values[nrow(hg19)]+hg19$length[nrow(hg19)]), ylim = c(0, 30), ylab = "Copy number", xlab = "", las=1, xaxt = "n", frame = FALSE, main = paste0("Normalized, GC-corrected, Filtered CN in ", w, "; # hq cancer cells = ", max(df_final$sort_rank)))
axis(side = 1, labels = F, at = hg19$add_values)
abline(v = hg19$add_values, col = "gray", lwd = 1/3)
tempdf$state[tempdf$state < 0] <- 0
tempdf$state[tempdf$state >= 11] <- 11
segments(x0 = tempdf$plot.start[tempdf$copy != "NaN" & tempdf$copy < 30], y0 = tempdf$copy[tempdf$copy != "NaN" & tempdf$copy < 30], x1 = tempdf$plot.end[tempdf$copy != "NaN" & tempdf$copy < 30], y1 = tempdf$copy[tempdf$copy != "NaN" & tempdf$copy < 30], col = modcol[tempdf$state[tempdf$copy != "NaN" & tempdf$copy < 30]+1], lwd = 2/3)
if (nrow(tempdf[tempdf$copy != "NaN" & tempdf$copy >= 30,]) != 0){
  segments(x0 = tempdf$plot.start[tempdf$copy != "NaN" & tempdf$copy >= 30], y0 = 30, x1 = tempdf$plot.end[tempdf$copy != "NaN" & tempdf$copy >= 30], y1 = 30, col = "navy", lwd = 3)
}
dev.off()

# Panel C, bottom
geneofi <- c("EGFR", "TERT", "CDK4")
celldf <- read.csv("./ordered_loadings.csv")
colnames(celldf) <- c("cell_id", "e4", "e3", "e1", "e2", "e5")

## Removal of noise
celldf$e1[celldf$e1 < 5] <- 0
celldf$e2[celldf$e2 < 5] <- 0
celldf$e3[celldf$e3 < 5] <- 0
celldf$e4[celldf$e4 < 5] <- 0
celldf$e5[celldf$e5 < 5] <- 0

for (i in c("e5", "e4", "e3", "e2", "e1")){
  celldf <- celldf[order(celldf[,i], decreasing = T),]
}
celldf <- celldf[celldf$cell_id %in% target_cell_ids[target_cell_ids %in% row.names(ad$obs)],]
celldf$sort_rank <- c(nrow(celldf):1)
celldf <- celldf[order(celldf$sort_rank, decreasing = T),]

df <- east[east$individual == w,]
setDT(df)
colnames(df)[2] <- "chr"
setDT(celldf)

df_final[celldf, sort_rank := i.sort_rank, on = "cell_id"]
df_final <- df_final[
  df,
  on = .(chr = chr,
         start >= start, 
         end <= end),
  
  j = .(
    chr = x.chr,
    start = x.start,
    end = x.end,
    reads = x.reads,
    copy = x.copy,
    state = x.state,
    cell_id = x.cell_id,
    normal_reads = x.normal_reads,
    bin_id = x.bin_id,
    segment_length = x.segment_length,
    plot.start = x.plot.start,
    plot.end = x.plot.end,
    centroacro = x.centroacro,
    blacklist = x.blacklist,
    xfilter = x.xfilter,
    sort_rank = x.sort_rank,
    regs = i.id
  ),
  
  nomatch = 0,
  allow.cartesian = TRUE
]

df_final[, chr_ord := {
  x <- sub("^chr", "", chr, ignore.case = TRUE)
  v <- suppressWarnings(as.integer(x))
  v[is.na(v) & toupper(x) == "X"] <- 23
  v[is.na(v) & toupper(x) == "Y"] <- 24
  v[is.na(v) & toupper(x) %in% c("M", "MT")] <- 25 
  v[is.na(v)] <- 9999 
  v
}]
setorder(df_final, sort_rank, chr_ord, start, end)
df_final[, chr_ord := NULL]

df[, chr_ord := {
  x <- sub("^chr", "", chr, ignore.case = TRUE)
  v <- suppressWarnings(as.integer(x))
  v[is.na(v) & toupper(x) == "X"] <- 23
  v[is.na(v) & toupper(x) == "Y"] <- 24
  v[is.na(v) & toupper(x) %in% c("M", "MT")] <- 25 
  v[is.na(v)] <- 9999 
  v
}]
setorder(df, chr_ord, start, end)
df[, chr_ord := NULL]
df$seglen <- df$end - df$start + 1
gaplen <- round(sum(df$seglen) * 0.05)
df$add_value <- 0
if (nrow(df) >= 2){
  for (i in 2:nrow(df)){
    df$add_value[i] <- df$add_value[i-1] + df$seglen[i-1] + gaplen
  }
}

regs <- df$id
for (i in 1:length(regs)){
  df_final$plot.start[df_final$regs == regs[i]] <- df_final$start[df_final$regs == regs[i]] - df$start[df$id == regs[i]] + df$add_value[df$id == regs[i]]
  df_final$plot.end[df_final$regs == regs[i]] <- df_final$end[df_final$regs == regs[i]] - df$start[df$id == regs[i]] + df$add_value[df$id == regs[i]]
}

pdf("Fig4c.bottom.pdf", width = 4, height = 4)
cell_number <- max(df_final$sort_rank)
plot(NULL, xlim = c(0, df$add_value[nrow(df)] + df$seglen[nrow(df)]), ylim = c(0, cell_number), ylab = "Cancer Cell Copy Number", xlab = "", las=1, xaxt = "n", frame = FALSE)
text(x = hg19$text_location[1:24], y = -max(df_final$sort_rank)*0.07, labels = hg19$chr[1:24], xpd=NA)
df_final$state[df_final$state < 0] <- 0
df_final$logcopy <- log10(df_final$state + 1)
df_final$linearcopy <- df_final$copy
df_final$linearcopy[df_final$linearcopy > 100] <- 100

rect(
  xleft = df_final$plot.start,
  ybottom = cell_number - df_final$sort_rank,
  xright = df_final$plot.end,
  ytop = cell_number - df_final$sort_rank + 1,
  border = NA,
  col = viridis::viridis_pal(option = "C", direction = 1)(100)[
    cut(
      x = df_final$linearcopy,
      breaks = seq(0, max(df_final$linearcopy, na.rm = TRUE), length.out = 101),
      labels = FALSE,
      include.lowest = TRUE
    )
  ]
)

for (i in geneofi){
  segments(x0 = genedf$start[gsub("-", ".", genedf$gene) == i] - df$start[grepl(i, df$oncogene)] + df$add_value[grepl(i, df$oncogene)], y0 = -0.05*k, x1 = genedf$end[gsub("-", ".", genedf$gene) == i] - df$start[grepl(i, df$oncogene)] + df$add_value[grepl(i, df$oncogene)], y1 = -0.05*k, lwd = 1, xpd = T, col = "darkblue")
  text(x = genedf$start[gsub("-", ".", genedf$gene) == i] - df$start[grepl(i, df$oncogene)] + df$add_value[grepl(i, df$oncogene)], y = -0.2*k, labels = i, xpd = T)
}

dev.off()

pdf("Fig4c.bottom.label.pdf", width = 4, height = 2)
image(
  z = t(matrix(1:100, nrow = 1)),
  col = viridis::viridis_pal(option = "C", direction = 1)(100),
  axes = FALSE,
  xlab = "", ylab = ""
)

axis(
  side = 1,
  at = c(0, 50, 100, 150, 200, 250)[c(0, 50, 100, 150, 200, 250) <= max(df_final$linearcopy, na.rm = TRUE)] / 
    max(df_final$linearcopy, na.rm = TRUE),
  labels = c(0, 50, 100, 150, 200, 250)[c(0, 50, 100, 150, 200, 250) <= max(df_final$linearcopy, na.rm = TRUE)],
  tick = TRUE,
  line = 0
)

mtext("Copy Number", side = 1, line = 2)
dev.off()

# Panel C, middle
mat <- as.matrix(celldf[, c(2:6)])
col_names <- colnames(mat) 
min_val <- min(mat)
max_val <- 100

col_list <- list(
  "e5" = circlize::colorRamp2(c(min_val, max_val), c("white", "#2378b5")), # Green
  "e1" = circlize::colorRamp2(c(min_val, max_val), c("white", "#f57f20")), # Red
  "e2" = circlize::colorRamp2(c(min_val, max_val), c("white", "#2fa148")), # Blue
  "e3" = circlize::colorRamp2(c(min_val, max_val), c("white", "#d52c28")),  # Orange
  "e4" = circlize::colorRamp2(c(min_val, max_val), c("white", "#9268ad"))  # Orange
)

h_list <- NULL

for (i in seq_along(col_names)) {
  nm <- col_names[i]
  
  h <- Heatmap(
    mat[, i, drop = FALSE],
    name = nm,
    col = col_list[[nm]],
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_order = 1:nrow(mat)
  )
  
  h_list <- h_list + h
}

pdf("Fig4c.middle.loading.pdf", width = 4, height = 4)
draw(h_list)
dev.off()

# Panel C, right
i <- "130115A-R23-C51"
i <- "130115A-R24-C19"
i <- "130115A-R17-C49"
z <- 4
regions <- "5-1000000-23100000;7-52000000-58000000;12-58000000-58500000"
geneofi <- c("EGFR", "TERT", "CDK4")

svsketch.dlp.sv.seg.cn.onecell.cssv(w, i, z, regions, geneofi)

# Panel D: conceptual illustration

# Function: clone-specific SV list
clone.specific.sv <- function(w, ...){
  cloneofi <- as.character(list(...))
  # Step 1
  if (w == "2765_2"){
    cndf <- read.csv("./results/sv_cna_clones_nov_14_v2/2765_2/tables/2765_2_clusters_0.7_1000.csv")
  } else if (w == "GBM0510_custom"){
    cndf <- read.csv("./metadata/custom_clones/GBM0510_custom.clones.csv")
  } else if (w == "NCI-H69_custom"){
    cndf <- read.csv("./metadata/custom_clones/NCI-H69_custom.clones.csv")
  } else if (w == "GBM0721_custom"){
    cndf <- read.csv("./metadata/custom_clones/GBM0721_custom.clones.csv")
  } else {
    cndf <- read.csv(paste0("./results/sv_gc_normalized_oct_18/", w, "/tables/", w, "_clusters_0.7_1000.csv"))
  }
  colnames(cndf)[2] <- "clone_id"
  cndf$cell_id <- str_sub(cndf$cell_id, -15)
  
  if (w == "GBM0510_custom"){
    qc <- fread(paste0("./metadata/doublet_invariant_filtered/GBM0510_qc.filtered.csv"))
    qc$cell_id <- str_sub(qc$cell_id, -15)
    qc <- qc[cell_id %in% cndf$cell_id]
    cndf <- cndf[cndf$cell_id %in% qc$cell_id,]
  } else if (w == "NCI-H69_custom"){
    qc <- fread(paste0("./metadata/doublet_invariant_filtered/NCI-H69_qc.filtered.csv"))
    qc$cell_id <- str_sub(qc$cell_id, -15)
    qc <- qc[cell_id %in% cndf$cell_id]
    cndf <- cndf[cndf$cell_id %in% qc$cell_id,]
  } else if (w == "GBM0721_custom"){
    qc <- fread(paste0("./metadata/doublet_invariant_filtered/GBM0721_qc.filtered.csv"))
    qc$cell_id <- str_sub(qc$cell_id, -15)
    qc <- qc[cell_id %in% cndf$cell_id]
    cndf <- cndf[cndf$cell_id %in% qc$cell_id,]
  } else {
    qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
    qc$cell_id <- str_sub(qc$cell_id, -15)
    cndf <- cndf[cndf$cell_id %in% qc$cell_id[qc$decision == "include"],]
  }
  
  # Step 2
  filelist <- list.files(
    path = file.path("./analysis/sc_sv_genotyping", w),
    pattern = "split_alt_counts.csv",
    full.names = TRUE,
    recursive = TRUE
  )
  
  svmat <- rbindlist(
    lapply(filelist, function(f) {
      dt <- fread(f)
      dt
    }),
    use.names = TRUE, # match by column names
    fill = TRUE # add missing columns with NA
  )
  svmat$cell_id <- str_sub(svmat$cell_id, -15)
  
  colnames(svmat) <- gsub("_right", "", gsub("_left", "", colnames(svmat)))
  if (w %in% c("GBM0510_custom", "NCI-H69_custom", "GBM0721_custom")){
    svmat <- svmat[cell_id %in% qc$cell_id]
  } else {
    svmat <- svmat[cell_id %in% qc$cell_id[qc$decision == "include"]]
  }
  mat <- as.matrix(svmat[, .SD, .SDcols = 2:ncol(svmat)])
  rownames(mat) <- svmat[["cell_id"]]
  
  filelist <- list.files(
    path = file.path("./analysis/sc_sv_genotyping", w),
    pattern = "span_alt_counts.csv",
    full.names = TRUE,
    recursive = TRUE
  )
  
  svmat <- rbindlist(
    lapply(filelist, function(f) {
      dt <- fread(f)
      dt
    }),
    use.names = TRUE, # match by column names
    fill = TRUE # add missing columns with NA
  )
  svmat$cell_id <- str_sub(svmat$cell_id, -15)
  
  colnames(svmat) <- gsub("_right", "", gsub("_left", "", colnames(svmat)))
  if (w %in% c("GBM0510_custom", "NCI-H69_custom", "GBM0721_custom")){
    svmat <- svmat[cell_id %in% qc$cell_id]
  } else {
    svmat <- svmat[cell_id %in% qc$cell_id[qc$decision == "include"]]
  }
  svmat <- svmat[, .SD, .SDcols = c("cell_id", names(svmat)[names(svmat) %in% colnames(mat)])]
  mat2 <- as.matrix(svmat[, .SD, .SDcols = 2:ncol(svmat)])
  rownames(mat2) <- svmat[["cell_id"]]
  
  # 1. Check row names
  if (!identical(rownames(mat), rownames(mat2))) {
    stop("Row names do not match between the two matrices.")
  }
  
  # 2. Check column names
  if (!identical(colnames(mat), colnames(mat2))) {
    stop("Column names do not match between the two matrices.")
  }
  
  # 3. Safe addition
  mat <- mat + mat2
  
  # 4. Load svdf
  if (w == "GBM0510_custom"){
    svdf <- fread(paste0("./analysis/hmf/gridss_somatic/GBM0510/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.annotated.bedpe"))
  } else if (w == "NCI-H69_custom"){
    svdf <- fread(paste0("./analysis/hmf/gridss_tumoronly/NCI-H69/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.annotated.bedpe"))
  } else if (w == "GBM0721_custom"){
    svdf <- fread(paste0("./analysis/hmf/gridss_somatic/GBM0721_custom/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.annotated.bedpe"))
  } else {
    svdf <- fread(paste0("./analysis/hmf/gridss_", ifelse(w %in% indi$individual[indi$normal_wgs == "yes"], "somatic/", "tumoronly/"), w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.annotated.bedpe"))
  }
  cndf <- cndf[cndf$cell_id %in% rownames(mat),]
  svdf[svname %in% colnames(mat[cndf$cell_id[cndf$clone_id %in% cloneofi], ])[colSums(mat[cndf$cell_id[cndf$clone_id %in% cloneofi], ]) > 1][!(colnames(mat[cndf$cell_id[cndf$clone_id %in% cloneofi], ])[colSums(mat[cndf$cell_id[cndf$clone_id %in% cloneofi], ]) > 1] %in% colnames(mat[cndf$cell_id[!(cndf$clone_id %in% cloneofi)], ])[colSums(mat[cndf$cell_id[!(cndf$clone_id %in% cloneofi)], ]) > 1])]]$svname
}



# Function: single-cell cnsv plotter, considering clone-specific SV
svsketch.dlp.sv.seg.cn.onecell.cssv <- function(w, i, z, regions, geneofi){
  options(scipen = 999)
  regs <- strsplit(regions, ";", fixed = T)[[1]]
  
  df <- as.data.frame(regs)
  df$chr <- ""
  df$start <- 0
  df$end <- 0
  for (j in 1:nrow(df)){
    df$chr[j] <- strsplit(df$regs[j], "-", fixed = T)[[1]][1]
    df$start[j] <- as.numeric(strsplit(df$regs[j], "-", fixed = T)[[1]][2])
    df$end[j] <- as.numeric(strsplit(df$regs[j], "-", fixed = T)[[1]][3])
  }
  df$seglen <- df$end - df$start + 1
  setDT(df)
  gaplen <- round(sum(df$seglen) * 0.05)
  df$add_value <- 0
  if (nrow(df) >= 2){
    for (j in 2:nrow(df)){
      df$add_value[j] <- df$add_value[j-1] + df$seglen[j-1] + gaplen
    }
  }
  
  df$gene <- ""
  for (j in geneofi){
    if (df$gene[df$chr == genedf$chromosome[genedf$gene == j] & df$start <= genedf$end[genedf$gene == j] & df$end >= genedf$start[genedf$gene == j]][1] != ""){
      df$gene[df$chr == genedf$chromosome[genedf$gene == j] & df$start <= genedf$end[genedf$gene == j] & df$end >= genedf$start[genedf$gene == j]][1] <- paste0(df$gene[df$start <= genedf$end[genedf$gene == j] & df$end >= genedf$start[genedf$gene == j]][1], ";", j)
    }
    df$gene[df$chr == genedf$chromosome[genedf$gene == j] & df$start <= genedf$end[genedf$gene == j] & df$end >= genedf$start[genedf$gene == j]][1] <- j
  }
  
  chrs <- df$chr
  
  if (w == "GBM0510_custom"){
    qc <- fread(paste0("./metadata/doublet_invariant_filtered/GBM0510_qc.filtered.csv"))
  } else if (w == "NCI-H69_custom"){
    qc <- fread(paste0("./metadata/doublet_invariant_filtered/NCI-H69_qc.filtered.csv"))
  } else if (w == "GBM0721_custom"){
    qc <- fread(paste0("./metadata/doublet_invariant_filtered/GBM0721_qc.filtered.csv"))
  } else {
    qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
  }
  qc$cell_id <- str_sub(qc$cell_id, -15)
  
  if (w == "GBM0510_custom"){
    df_final <- fread(paste0("./analysis/normalization/cn_dataframe/GBM0510.normalized.gc_corrected.filters.annotated.csv"))
  } else if (w == "NCI-H69_custom"){
    df_final <- fread(paste0("./analysis/normalization/cn_dataframe/NCI-H69.normalized.gc_corrected.filters.annotated.csv"))
  } else if (w == "GBM0721_custom"){
    df_final <- fread(paste0("./analysis/normalization/cn_dataframe/GBM0721.normalized.gc_corrected.filters.annotated.csv"))
  } else {
    df_final <- fread(paste0("./analysis/normalization/cn_dataframe/", w, ".normalized.gc_corrected.filters.annotated.csv"))
  }
  df_final$cell_id <- str_sub(df_final$cell_id, -15)
  
  if (w == "GBM0510_custom"){
    celldf <- read.csv("./analysis/metadata/custom_clones/GBM0510_custom.clones.csv")
    df_final <- df_final[
      normal_reads != 0 &
        segment_length >= 1000 &
        blacklist == "no" &
        xfilter == "no" &
        cell_id %in% celldf$cell_id
    ]
  } else if (w == "NCI-H69_custom"){
    celldf <- read.csv("./analysis/metadata/custom_clones/NCI-H69_custom.clones.csv")
    df_final <- df_final[
      normal_reads != 0 &
        segment_length >= 1000 &
        blacklist == "no" &
        xfilter == "no" &
        cell_id %in% celldf$cell_id
    ]
  } else if (w == "GBM39-DM"){
    celldf <- read.csv("./analysis/metadata/custom_clones/GBM39-DM_custom.clones.csv")
    df_final <- df_final[
      normal_reads != 0 &
        segment_length >= 1000 &
        blacklist == "no" &
        xfilter == "no" &
        cell_id %in% celldf$cell_id
    ]
  } else if (w == "GBM0721_custom"){
    celldf <- read.csv("./analysis/metadata/custom_clones/GBM0721_custom.clones.csv")
    df_final <- df_final[
      normal_reads != 0 &
        segment_length >= 1000 &
        blacklist == "no" &
        xfilter == "no" &
        cell_id %in% celldf$cell_id
    ]
  } else {
    df_final <- df_final[
      normal_reads != 0 &
        segment_length >= 1000 &
        blacklist == "no" &
        xfilter == "no" &
        cell_id %in% qc$cell_id[qc$decision == "include"]
    ]
  }
  setDT(df_final)
  setDT(df)
  df_final[, chr := as.character(chr)]
  df[, chr := as.character(chr)]
  df_final[, .rowid := .I]
  
  setkey(df_final, chr, start, end)
  setkey(df, chr, start, end)
  idx <- foverlaps(df_final, df, type = "any", which = TRUE, nomatch = 0L)
  df_final <- df_final[unique(idx$xid)]
  
  df_final_trimmed <- df_final[
    df,
    on = .(chr = chr,
           start <= end, 
           end >= start), 
    
    j = .(
      chr = x.chr,
      start = pmax(x.start, i.start), 
      end = pmin(x.end, i.end),
      regs = i.regs,
      reads = x.reads,
      copy = x.copy,
      state = x.state,
      cell_id = x.cell_id,
      normal_reads = x.normal_reads,
      bin_id = x.bin_id,
      segment_length = x.segment_length,
      plot.start = x.plot.start,
      plot.end = x.plot.end,
      centroacro = x.centroacro,
      blacklist = x.blacklist,
      xfilter = x.xfilter,
      sort_rank = x.sort_rank,
      .rowid = x..rowid 
    ),
    
    nomatch = 0, 
    allow.cartesian = TRUE 
  ]
  
  df_final <- df_final_trimmed
  
  df_final[, chr_ord := {
    x <- sub("^chr", "", chr, ignore.case = TRUE)
    v <- suppressWarnings(as.integer(x))
    v[is.na(v) & toupper(x) == "X"] <- 23
    v[is.na(v) & toupper(x) == "Y"] <- 24
    v[is.na(v) & toupper(x) %in% c("M", "MT")] <- 25 # optional
    v[is.na(v)] <- 9999 # push any non-standard names to the end
    v
  }]
  setorder(df_final, chr_ord, start, end)
  df_final[, chr_ord := NULL]
  df_final[df, on = .(regs),
           `:=`(
             start = pmax(start, i.start),
             end   = pmin(end,   i.end)
           )]
  
  
  df[, chr_ord := {
    x <- sub("^chr", "", chr, ignore.case = TRUE)
    v <- suppressWarnings(as.integer(x))
    v[is.na(v) & toupper(x) == "X"] <- 23
    v[is.na(v) & toupper(x) == "Y"] <- 24
    v[is.na(v) & toupper(x) %in% c("M", "MT")] <- 25 # optional
    v[is.na(v)] <- 9999 # push any non-standard names to the end
    v
  }]
  setorder(df, chr_ord, start, end)
  df[, chr_ord := NULL]
  
  if (i %in% unique(df_final$cell_id)){
    for (j in 1:length(regs)){
      df_final$plot.start[df_final$regs == regs[j]] <- df_final$start[df_final$regs == regs[j]] - df$start[df$regs == regs[j]] + df$add_value[df$regs == regs[j]]
      df_final$plot.end[df_final$regs == regs[j]] <- df_final$end[df_final$regs == regs[j]] - df$start[df$regs == regs[j]] + df$add_value[df$regs == regs[j]]
    }
    
    if (w == "GBM0510_custom"){
      svdf <- fread(paste0("./analysis/hmf/gridss_somatic/GBM0510/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.annotated.bedpe"))
    } else if (w == "NCI-H69_custom"){
      svdf <- fread(paste0("./analysis/hmf/gridss_tumoronly/NCI-H69/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.annotated.bedpe"))
    } else if (w == "GBM0721_custom"){
      svdf <- fread(paste0("./analysis/hmf/gridss_somatic/GBM0721_custom/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.annotated.bedpe"))
    } else {
      svdf <- fread(paste0("./analysis/hmf/gridss_", ifelse(w %in% indi$individual[indi$normal_wgs == "yes"], "somatic/", "tumoronly/"), w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.annotated.bedpe"))
    }
    svdf[, chr1 := as.character(chr1)]
    svdf[, chr2 := as.character(chr2)]
    swap <- chr_rank(svdf$chr2) < chr_rank(svdf$chr1)
    svdf[swap, c("chr1","start1","end1","ori1","chr2","start2","end2","ori2") :=
           .(chr2, start2, end2, ori2, chr1, start1, end1, ori1)]
    
    svdf$pos1 <- round((svdf$start1 + svdf$end1)/2)
    svdf$pos2 <- round((svdf$start2 + svdf$end2)/2)
    colnames(svdf)[colnames(svdf) == "ori1"] <- "strand1"
    colnames(svdf)[colnames(svdf) == "ori2"] <- "strand2"
    
    i1 <- svdf[df, on=.(chr1=chr, start1<=end, end1>=start), which=TRUE, nomatch=0L]
    i2 <- svdf[df, on=.(chr2=chr, start2<=end, end2>=start), which=TRUE, nomatch=0L]
    svdf <- svdf[unique(c(i1, i2))]
    svdf[, chr_ord := {
      x <- sub("^chr", "", chr1, ignore.case = TRUE)
      v <- suppressWarnings(as.integer(x))
      v[is.na(v) & toupper(x) == "X"] <- 23
      v[is.na(v) & toupper(x) == "Y"] <- 24
      v[is.na(v) & toupper(x) %in% c("M", "MT")] <- 25 # optional
      v[is.na(v)] <- 9999 # push any non-standard names to the end
      v
    }]
    setorder(svdf, chr_ord, start1, end1)
    svdf[, chr_ord := NULL]
    
    svdf[, reg1 := df[
      # Use .SD for all rows of svdf
      .SD, 
      # Non-equi join conditions
      on = .(chr = chr1, start <= pos1, end >= pos1), 
      j = regs, # Retrieve the 'regs' column from df
      nomatch = NA_character_
      # REMOVE THE REDUNDANT ]$regs HERE
    ]]
    
    # Similarly for the second join:
    svdf[, reg2 := df[
      .SD, 
      on = .(chr = chr2, start <= pos2, end >= pos2), 
      j = regs, 
      nomatch = NA_character_
    ]]
    
    filelist <- list.files(
      path = file.path("./analysis/sc_sv_genotyping", w),
      pattern = "split_alt_counts.csv",
      full.names = TRUE,
      recursive = TRUE
    )
    
    svmat <- rbindlist(
      lapply(filelist, function(f) {
        dt <- fread(f)
        dt
      }),
      use.names = TRUE, 
      fill = TRUE 
    )
    svmat$cell_id <- str_sub(svmat$cell_id, -15)
    
    colnames(svmat) <- gsub("_right", "", gsub("_left", "", colnames(svmat)))
    if (w %in% c("GBM0510_custom", "NCI-H69_custom")){
      svmat <- svmat[cell_id %in% celldf$cell_id]
    } else {
      svmat <- svmat[cell_id %in% qc$cell_id[qc$decision == "include"]]
    }
    mat <- as.matrix(svmat[, .SD, .SDcols = 2:ncol(svmat)])
    rownames(mat) <- svmat[["cell_id"]]
    
    filelist <- list.files(
      path = file.path("./analysis/sc_sv_genotyping", w),
      pattern = "span_alt_counts.csv",
      full.names = TRUE,
      recursive = TRUE
    )
    
    svmat <- rbindlist(
      lapply(filelist, function(f) {
        dt <- fread(f)
        dt
      }),
      use.names = TRUE, 
      fill = TRUE 
    )
    svmat$cell_id <- str_sub(svmat$cell_id, -15)
    
    colnames(svmat) <- gsub("_right", "", gsub("_left", "", colnames(svmat)))
    if (w %in% c("GBM0510_custom", "NCI-H69_custom", "GBM0721_custom")){
      svmat <- svmat[cell_id %in% celldf$cell_id]
    } else {
      svmat <- svmat[cell_id %in% qc$cell_id[qc$decision == "include"]]
    }
    svmat <- svmat[, .SD, .SDcols = c("cell_id", names(svmat)[names(svmat) %in% colnames(mat)])]
    mat2 <- as.matrix(svmat[, .SD, .SDcols = 2:ncol(svmat)])
    rownames(mat2) <- svmat[["cell_id"]]
    
    if (!identical(rownames(mat), rownames(mat2))) {
      stop("Row names do not match between the two matrices.")
    }
    
    if (!identical(colnames(mat), colnames(mat2))) {
      stop("Column names do not match between the two matrices.")
    }
    
    mat <- mat + mat2
    svdf <- svdf[svname %in% colnames(mat)]
    
    if (w %in% c("GBM0510_custom", "NCI-H69_custom", "GBM0721_custom")){
      celllist <- celldf$cell_id
    } else {
      celllist <- qc$cell_id[qc$decision == "include" & qc$cell_id %in% rownames(mat)]
    }
    
    pdf(paste0("figures/svsketch/sv-seg-cn/", w, ".", i, ".", paste(chrs, collapse = "-"), ".dlp.sv-seg-cn.onecell.cssv.svsketch.pdf"), width=z, height=3.7) 
    cnvdf <- as.data.frame(df_final[cell_id == i])
    cnvdf <- cnvdf[,c("chr", "plot.start", "plot.end", "copy", "state")]
    k <- max(cnvdf$copy)
    colnames(cnvdf) <- c("chromosome", "plotstart", "plotend", "CopyNumber", "state")
    cnvdf$state[cnvdf$state >= 11] <- 11
    cnvdf$state[cnvdf$state < 0] <- 0
    cnvdf$plotstart <- as.numeric(cnvdf$plotstart)
    cnvdf$plotend <- as.numeric(cnvdf$plotend)
    
    ## PLOTTING COPY NUMBERS
    if (k < 5){
      k <- 5
    }
    plot(0, type='n', xlim=c(0, df$add_value[nrow(df)] + df$seglen[nrow(df)]), ylim=c(0,1.1*k), xaxt = 'n', frame=F, las=1, xlab=paste0("Position on chromosome ", paste(unique(df$chr), collapse = "-"), " in Mbp"), ylab="Total copy number", main = paste0(i, " quality = ", qc$quality[qc$cell_id == i]))
    for (j in 1:nrow(df)){
      axis(1, at = c(df$add_value[j], df$add_value[j] + df$seglen[j]), tick = c(df$add_value[j]+1, df$add_value[j] + df$seglen[j]), labels = c(df$start[j], df$end[j]), las=2)
    }
    segments(cnvdf$plotstart, cnvdf$CopyNumber, cnvdf$plotend, cnvdf$CopyNumber, col=modcol[as.numeric(cnvdf$state)+1], lwd=2/3)
    
    # GENES OF INTEREST
    for (j in geneofi){
      segments(x0 = genedf$start[genedf$gene == j] - df$start[grepl(j, df$gene)] + df$add_value[grepl(j, df$gene)], y0 = 1.1*k, x1 = genedf$end[genedf$gene == j] - df$start[grepl(j, df$gene)] + df$add_value[grepl(j, df$gene)], y1 = 1.1*k, col = "darkblue", lwd = 2, xpd = T)
      text(x = genedf$start[genedf$gene == j] - df$start[grepl(j, df$gene)] + df$add_value[grepl(j, df$gene)], y = 1.2*k, labels = j, xpd = T)
    }
    
    l <- 5
    
    if (length(colnames(mat)[mat[i,] != 0]) != 0 & nrow(svdf[svname %in% colnames(mat)[mat[i,] != 0] & !(type == "DEL" & svlen < 500)]) != 0){
      isv <- as.data.frame(svdf[svname %in% colnames(mat)[mat[i,] != 0] & !(type == "DEL" & svlen < 500)])
      isv$pos1 <- as.numeric(isv$pos1)
      isv$pos2 <- as.numeric(isv$pos2)
      isv$plotpos1 <- isv$pos1
      isv$plotpos2 <- isv$pos2
      
      for (j in 1:nrow(df)){
        isv$plotpos1[isv$reg1 == df$regs[j] & !is.na(isv$reg1)] <- isv$plotpos1[isv$reg1 == df$regs[j] & !is.na(isv$reg1)] - df$start[j] + df$add_value[j]
        isv$plotpos2[isv$reg2 == df$regs[j] & !is.na(isv$reg2)] <- isv$plotpos2[isv$reg2 == df$regs[j] & !is.na(isv$reg2)] - df$start[j] + df$add_value[j]
      }
      
      if (max(mat[i,svdf[svname %in% colnames(mat)[mat[i,] != 0] & !(type == "DEL" & svlen < 500)]$svname]) > 5){
        l <- max(mat[i,svdf[svname %in% colnames(mat)[mat[i,] != 0] & !(type == "DEL" & svlen < 500)]$svname])
      }
      
      if (w == "2765_2"){
        cndf <- read.csv("./results/sv_cna_clones_nov_14_v2/2765_2/tables/2765_2_clusters_0.7_1000.csv")
      } else if (w == "GBM0510_custom"){
        cndf <- read.csv("./analysis/metadata/custom_clones/GBM0510_custom.clones.csv")
      } else if (w == "NCI-H69_custom"){
        cndf <- read.csv("./analysis/metadata/custom_clones/NCI-H69_custom.clones.csv")
      } else if (w == "GBM39-DM"){
        cndf <- read.csv(paste0("./results/sv_gc_normalized_oct_18/", w, "/tables/", w, "_clusters_0.7_1000.csv"))
      }
      colnames(cndf)[2] <- "clone_id"
      cndf$cell_id <- str_sub(cndf$cell_id, -15)
      
      isv$clonespecific <- "no"
      if (w == "NCI-H69_custom"){
        if (i %in% cndf$cell_id){
          isv$clonespecific[isv$svname %in% clone.specific.sv(w, "A", cndf$clone_id[cndf$cell_id == i])] <- "yes"
        }
      } else {
        if (i %in% cndf$cell_id){
          isv$clonespecific[isv$svname %in% clone.specific.sv(w, cndf$clone_id[cndf$cell_id == i])] <- "yes"
        }
      }
      
      if (nrow(isv[isv$chr1 %in% chrs & isv$chr2 %in% chrs,]) != 0){
        isvin <- isv[isv$chr1 %in% chrs & isv$chr2 %in% chrs,]
        par(new = T)
        plot(0, type='n', xlim=c(0, df$add_value[nrow(df)] + df$seglen[nrow(df)]), ylim=c(0,1.1*l), frame=F, axes = F, xlab = "", ylab = "")
        axis(4, las=1)
        mtext(paste0("SV support in cell ", i), side=4, line=3)
        
        theta=seq(0,pi, len=100)
        isvin$rad=abs(isvin$plotpos2-isvin$plotpos1)/2
        
        for (j in 1:nrow(isvin)){
          if (isvin$strand1[j] == "+" & isvin$strand2[j] == "-"){
            svcol <- ifelse(isvin$clonespecific[j] == "yes", rgb(0/255,0/255,255/255,1), rgb(0/255,0/255,255/255,.3))
            x <- isvin$rad[j]*cos(theta)+((isvin$plotpos1[j]+isvin$plotpos2[j])/2)
            y <- l*0.1*sin(theta)+mat[i,isvin$svname[j]]
          } else if (isvin$strand1[j] == "-" & isvin$strand2[j] == "+"){
            svcol <- ifelse(isvin$clonespecific[j] == "yes", rgb(0/255,128/255,128/255,1), rgb(0/255,128/255,128/255,.3))
            x <- isvin$rad[j]*cos(theta)+((isvin$plotpos1[j]+isvin$plotpos2[j])/2)
            y <- -l*0.1*sin(theta)+mat[i,isvin$svname[j]]
          } else if (isvin$strand1[j] == "-" & isvin$strand2[j] == "-"){
            svcol <- ifelse(isvin$clonespecific[j] == "yes", rgb(220/255,20/255,60/255,1), rgb(220/255,20/255,60/255,.3))
            x <- isvin$rad[j]*cos(theta)+((isvin$plotpos1[j]+isvin$plotpos2[j])/2)
            y <- l*0.1*sin(theta)+mat[i,isvin$svname[j]]
          } else if (isvin$strand1[j] == "+" & isvin$strand2[j] == "+"){
            svcol <- ifelse(isvin$clonespecific[j] == "yes", rgb(128/255,128/255,0/255,1), rgb(128/255,128/255,0/255,.3))
            x <- isvin$rad[j]*cos(theta)+((isvin$plotpos1[j]+isvin$plotpos2[j])/2)
            y <- -l*0.1*sin(theta)+mat[i,isvin$svname[j]]
          }
          lines(x,y,col=svcol, xpd=TRUE, lwd=2/3)
          segments(isvin$plotpos1[j], 0, isvin$plotpos1[j], mat[i,isvin$svname[j]], lty=1, xpd=TRUE, col=svcol, lwd=2/3)
          segments(isvin$plotpos2[j], 0, isvin$plotpos2[j], mat[i,isvin$svname[j]], lty=1, xpd=TRUE, col=svcol, lwd=2/3)
        }
      }
      
      if (nrow(isv[!(isv$chr1 %in% chrs) | !(isv$chr2 %in% chrs),]) != 0){
        isvout <- isv[!(isv$chr1 %in% chrs) | !(isv$chr2 %in% chrs),]
        for (j in 1:nrow(isvout)){
          svcol <- ifelse(isvout$clonespecific[j] == "yes", rgb(85/255,26/255,139/255,1), rgb(85/255,26/255,139/255,.15))
          if (isvout$chr1[j] %in% chrs){
            arrows(isvout$plotpos1[j], 0, isvout$plotpos1[j], mat[i,isvout$svname[j]], lty=1, xpd=TRUE, col=svcol, code=2, length = 0.1, lwd=2/3)
          } else {
            arrows(isvout$plotpos2[j], 0, isvout$plotpos2[j], mat[i,isvout$svname[j]], lty=1, xpd=TRUE, col=svcol, code=2, length = 0.1, lwd=2/3)
          }
        }
      }
    }
    dev.off()
  }
}


# Function: pseudobulk cnsv plotter
svsketch.dlp.pseudobulk <- function(i, chrs, geneofi){
  chrs <- strsplit(chrs, ";", fixed = T)[[1]]
  geneofi <- strsplit(geneofi, ";", fixed = T)[[1]]
  
  
  qc <- fread(paste0("./metadata/doublet_invariant_filtered/", i, "_qc.filtered.csv"))
  qc$cell_id <- str_sub(qc$cell_id, -15)
  
  df_final <- fread(paste0("./normalization/cn_dataframe/", i, ".normalized.gc_corrected.filters.annotated.csv"))
  df_final$cell_id <- str_sub(df_final$cell_id, -15)
  
  df_final <- df_final[
    normal_reads != 0 &
      segment_length >= 1000 &
      blacklist == "no" &
      xfilter == "no" &
      cell_id %in% qc$cell_id[qc$decision == "include"]
  ]
  
  setDT(df_final)
  df_final <- df_final[chr %in% chrs]
  cnvdf <- consensuscopynumber(as.data.frame(df_final[,c("chr", "start", "end", "reads", "copy", "state", "cell_id")]))
  
  gaplen <- sum(hg19$length[hg19$chr %in% chrs]) * 0.05
  deltas <- c()
  dv <- 0
  deltas[1] <- dv
  if (length(chrs) > 1){
    for (j in 2:length(chrs)){
      dv <- dv + hg19$length[hg19$chr == chrs[j-1]] + gaplen
      deltas[j] <- dv
    }
  }
  deltadf <- as.data.frame(chrs)
  deltadf$delta <- deltas
  deltadf$chrs <- as.character(deltadf$chrs)
  
  cnvdf <- cnvdf[cnvdf$chr %in% chrs,]
  for (j in 1:nrow(cnvdf)){
    cnvdf$start[j] <- cnvdf$start[j] + deltadf$delta[deltadf$chrs == cnvdf$chr[j]]
    cnvdf$end[j] <- cnvdf$end[j] + deltadf$delta[deltadf$chrs == cnvdf$chr[j]]
  }
  cnvdf <- cnvdf[order(cnvdf$start, decreasing = F),]
  cnvdf$state[cnvdf$state > 11] <- 11
  
  
  isv <- read.csv(paste0("analysis/hmf/gridss_", ifelse(indi$normal_wgs[indi$individual == i] == "yes", "somatic", "tumoronly"), "/", i, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.annotated.bedpe"), header=T, as.is=T, sep="\t")
  isv$pos1 <- round((isv$start1 + isv$end1)/2)
  isv$pos2 <- round((isv$start2 + isv$end2)/2)
  isv$svlen <- NA
  isv$svlen[isv$chr1 == isv$chr2] <- isv$pos2[isv$chr1 == isv$chr2] - isv$pos1[isv$chr1 == isv$chr2]
  
  isv <- isv[isv$chr1 %in% chrs | isv$chr2 %in% chrs,]
  isv <- isv[!is.na(isv$vf),]
  isv <- isv[!(isv$chr1 == isv$chr2 & isv$type == "DEL" & isv$svlen <= 100),]
  isv <- isv[!(isv$type == "TRA" & isv$vf <= 4),]
  isv <- isv[!(isv$type %in% c("t2tFBI", "h2hFBI") & isv$vf <= 4),]
  colnames(isv)[colnames(isv) == "ori1"] <- "strand1"
  colnames(isv)[colnames(isv) == "ori2"] <- "strand2"
  
  if (nrow(isv) == 0){
    return(NULL)
  }
  
  for (j in 1:nrow(isv)){
    if (isv$chr1[j] %in% chrs){
      isv$pos1[j] <- isv$pos1[j] + deltadf$delta[deltadf$chrs == isv$chr1[j]]
    }
    if (isv$chr2[j] %in% chrs){
      isv$pos2[j] <- isv$pos2[j] + deltadf$delta[deltadf$chrs == isv$chr2[j]]
    }
  }
  
  theta=seq(0,pi, len=100)
  isvin <- isv[isv$chr1 %in% chrs & isv$chr2 %in% chrs,]
  isvin$rad=abs(isvin$pos2-isvin$pos1)/2
  isvout <- isv[!(isv$chr1 %in% chrs) | !(isv$chr2 %in% chrs),]
  max_val <- max(isv$vf)
  
  pdf(paste0("Fig4a.", i, "_", paste(geneofi, collapse = "-"), ".pdf"), width=ifelse(length(chrs) == 1, 4, 8), height=4)
  par(mar=c(5.1,4.1,4.1,4.1)) # default
  
  plot(0, type='n', xlim=c(0, hg19$length[hg19$chr == deltadf$chrs[nrow(deltadf)]] + deltadf$delta[nrow(deltadf)]), ylim=c(0,1.1*max_val), frame=F, axes = F, xlab = "", ylab = "")
  axis(4, las=1)
  
  if (nrow(isvout) != 0){
    for (j in 1:nrow(isvout)){
      svcol <- rgb(85/255,26/255,139/255,.15)
      if (isvout$chr1[j] %in% chrs){
        arrows(isvout$pos1[j], 0, isvout$pos1[j], isvout$vf[j], lty=1, xpd=TRUE, col=svcol, code=2, lwd=2/3, length=0.05)
      } 
      if (isvout$chr2[j] %in% chrs) {
        arrows(isvout$pos2[j], 0, isvout$pos2[j], isvout$vf[j], lty=1, xpd=TRUE, col=svcol, code=2, lwd=2/3, length=0.05)
      }
    }
  }
  
  for (j in 1:nrow(isvin)){
    if (isvin$strand1[j] == "+" & isvin$strand2[j] == "-"){
      svcol <- rgb(0/255,0/255,255/255,.3)
      x <- isvin$rad[j]*cos(theta)+((isvin$pos1[j]+isvin$pos2[j])/2)
      y <- max_val*0.1*sin(theta)+isvin$vf[j]
    } else if (isvin$strand1[j] == "-" & isvin$strand2[j] == "+"){
      svcol <- rgb(0/255,128/255,128/255,.3)
      x <- isvin$rad[j]*cos(theta)+((isvin$pos1[j]+isvin$pos2[j])/2)
      y <- -max_val*0.1*sin(theta)+isvin$vf[j]
    } else if (isvin$strand1[j] == "+" & isvin$strand2[j] == "+"){
      if (isvin$type[j] == "t2tFBI"){
        svcol <- rgb(220/255,20/255,60/255,.6)
        x <- 0.025*sum(hg19$length[hg19$chr %in% chrs])*cos(seq(0,1.9*pi, len=100))+((isvin$pos1[j]+isvin$pos2[j])/2)-0.025*sum(hg19$length[hg19$chr %in% chrs])
        y <- max_val*0.05*sin(seq(0,1.9*pi, len=100))+isvin$vf[j]
      } else {
        svcol <- rgb(220/255,20/255,60/255,.3)
        x <- isvin$rad[j]*cos(theta)+((isvin$pos1[j]+isvin$pos2[j])/2)
        y <- max_val*0.1*sin(theta)+isvin$vf[j]
      }
    } else if (isvin$strand1[j] == "-" & isvin$strand2[j] == "-"){
      if (isvin$type[j] == "h2hFBI"){
        svcol <- rgb(128/255,128/255,0/255,.6)
        x <- 0.025*sum(hg19$length[hg19$chr %in% chrs])*-cos(seq(0,1.9*pi, len=100))+((isvin$pos1[j]+isvin$pos2[j])/2)+0.025*sum(hg19$length[hg19$chr %in% chrs])
        y <- max_val*0.05*sin(seq(0,1.9*pi, len=100))+isvin$vf[j]
      } else {
        svcol <- rgb(128/255,128/255,0/255,.3)
        x <- isvin$rad[j]*cos(theta)+((isvin$pos1[j]+isvin$pos2[j])/2)
        y <- -max_val*0.1*sin(theta)+isvin$vf[j]
      }
    }
    lines(x,y,col=svcol, xpd=TRUE, lwd=2/3)
    segments(isvin$pos1[j], 0, isvin$pos1[j], isvin$vf[j], lty=1, xpd=TRUE, col=svcol, lwd=2/3)
    segments(isvin$pos2[j], 0, isvin$pos2[j], isvin$vf[j], lty=1, xpd=TRUE, col=svcol, lwd=2/3)
  }
  
  ## Plotting copy numbers
  par(new = T)
  
  k <- max(cnvdf$copy)
  if (k < 5){
    k <- 5
  }
  plot(0, type='n', xlim=c(0, max(cnvdf$end)), ylim=c(0,1.1*k), frame=F, xaxt="n", las=1, xlab=paste0("Position on chromosomes ", paste(chrs, collapse = ", ")), ylab="Allelic copy number")
  
  tickloc <- c()
  ticklab <- c()
  numlist <- NULL
  lablist <- NULL
  for (j in 1:nrow(deltadf)){
    numlist <- seq(0, floor(hg19$length[hg19$chr == deltadf$chrs[j]]/20000000)*20000000, by = 20000000) + deltadf$delta[j]
    tickloc <- c(tickloc, numlist)
    lablist <- seq(0, floor(hg19$length[hg19$chr == deltadf$chrs[j]]/20000000)*20, by = 20)
    ticklab <- c(ticklab, lablist)
  }
  axis(1, at = tickloc, labels = ticklab)
  
  segments(cnvdf$start, cnvdf$copy, cnvdf$end, cnvdf$copy, col = modcol[cnvdf$state+1], lwd = 4/3, lend = 2)
  
  ## Genes of interest
  for (j in geneofi){
    if (genedf$chromosome[gsub("-", ".", genedf$gene) == j] %in% chrs){
      segments(x0 = genedf$start[gsub("-", ".", genedf$gene) == j] + deltadf$delta[deltadf$chrs == genedf$chromosome[gsub("-", ".", genedf$gene) == j]], y0 = 1.1*k, x1 = genedf$end[gsub("-", ".", genedf$gene) == j] + deltadf$delta[deltadf$chrs == genedf$chromosome[gsub("-", ".", genedf$gene) == j]], y1 = 1.1*k, col = "darkgreen", lwd = 2, xpd = T)
      text(x = genedf$start[gsub("-", ".", genedf$gene) == j] + deltadf$delta[deltadf$chrs == genedf$chromosome[gsub("-", ".", genedf$gene) == j]], y = 1.2*k, labels = j, xpd = T)
    }
  }
  
  ## Ideogram
  for (j in chrs){
    rect(cytoband$start[cytoband$chr == j] + deltadf$delta[deltadf$chrs == j], (-0.22)*k, cytoband$end[cytoband$chr == j] + deltadf$delta[deltadf$chrs == j], (-0.17)*k, col = cytoband$color[cytoband$chr == j], xpd = T, lwd = 0.5)
  }
  
  dev.off()
}