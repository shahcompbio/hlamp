# Figure 3
library(stringr)
library(ComplexHeatmap)
library(circlize)
options(scipen = 999)
indi <- read.csv("./metadata/dlp.cohort.summary.individual_111925.tsv", header=T, as.is=T, sep="\t")
east <- read.csv("./modeling/amplified_regions_summary_table_011026.csv")

# Panel A, left
w <- "GBM0721"

pdf(paste0("Fig3a.left.", w, ".pdf"), height=4, width=4.5)
layout(matrix(1:2, nrow=2), heights=c(0.4, 0.6))

## 1. Top pseudobulk CN plot
qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)
tempdf <- read.csv(paste0("./normalization/subclonal_cn_v4/", w, "_clusters_0.7_1000_subclonal_cn.gc_corrected.csv"))

par(mar = c(0.5, 4.1, 4.1, 2.1))
plot(NULL, xlim = c(hg19$add_values[1], hg19$add_values[nrow(hg19)]+hg19$length[nrow(hg19)]), ylim = c(0, 30), ylab = "Copy number", xlab = "", las=1, xaxt = "n", frame = FALSE, main = paste0(w, "; # hq cancer cells = ", indi$final_cancer_cell_ct[indi$individual == w]))
axis(side = 1, labels = F, at = hg19$add_values)
#text(x = hg19$text_location[1:24], y = -3, labels = hg19$chr[1:24], xpd=NA)
abline(v = hg19$add_values, col = "gray")
tempdf$state[tempdf$state < 0] <- 0
tempdf$state[tempdf$state >= 11] <- 11
segments(x0 = tempdf$plot.start[tempdf$copy != "NaN" & tempdf$copy < 30], y0 = tempdf$copy[tempdf$copy != "NaN" & tempdf$copy < 30], x1 = tempdf$plot.end[tempdf$copy != "NaN" & tempdf$copy < 30], y1 = tempdf$copy[tempdf$copy != "NaN" & tempdf$copy < 30], col = modcol[tempdf$state[tempdf$copy != "NaN" & tempdf$copy < 30]+1], lwd = 1)
if (nrow(tempdf[tempdf$copy != "NaN" & tempdf$copy >= 30,]) != 0){
  segments(x0 = tempdf$plot.start[tempdf$copy != "NaN" & tempdf$copy >= 30], y0 = 30, x1 = tempdf$plot.end[tempdf$copy != "NaN" & tempdf$copy >= 30], y1 = 30, col = "navy", lwd = 1)
}

## 2. Bottom single-cell CN heatmap 
df_final <- fread(paste0("./normalization/cn_dataframe/", w, ".normalized.gc_corrected.filters.annotated.csv"))
df_final$cell_id <- str_sub(df_final$cell_id, -15)

df_final <- df_final[
  normal_reads != 0 &
    segment_length >= 1000 &
    centroacro == "no" &
    blacklist == "no" &
    xfilter == "no" &
    cell_id %in% qc$cell_id[qc$decision == "include"]
]

par(mar = c(3, 4.1, 0.1, 2.1))
cell_number <- max(df_final$sort_rank)
plot(NULL, xlim = c(hg19$add_values[1], hg19$add_values[nrow(hg19)]+hg19$length[nrow(hg19)]), ylim = c(0, cell_number), ylab = "Cancer Cell Copy Number", xlab = "", las=1, xaxt = "n", frame = FALSE)
axis(side = 1, labels = F, at = hg19$add_values)
text(x = hg19$text_location[1:24], y = -max(df_final$sort_rank)*0.07, labels = hg19$chr[1:24], xpd=NA)
abline(v = hg19$add_values, col = "gray")
df_final$state[df_final$state < 0] <- 0
df_final$state[df_final$state >= 11] <- 11
rect(df_final$plot.start, cell_number-df_final$sort_rank, df_final$plot.end, cell_number-df_final$sort_rank+1, border = NA, col = modcol[df_final$state+1])
rm(df_final)
dev.off()

# Panel A, right
w <- "GBM0721"
i <- "EGFR" # x axis
j <- "MDM4" # y axis
coli <- rgb(39/255, 64/255, 139/255, .3)

pdf(paste0("fig3a.right.", w, ".", i, ".", j, ".pdf"), height=4, width=4)
cndf <- read.csv(paste0("./normalization/oncogene_cn_final/", w, ".sv-seg-cn.oncogene.cn.csv"))
cndf$cell_id <- str_sub(cndf$cell_id, -15)
metrics <- read.csv(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
metrics$cell_id <- str_sub(metrics$cell_id, -15)
cndf <- cndf[cndf$cell_id %in% metrics$cell_id[metrics$decision == "include"],]

max_val <- (floor(max(c(max(cndf[,j], cndf[,i])))/50)+1)*50
layout(matrix(c(2,4,1,3), 2, 2, byrow = TRUE), widths = c(4,1), heights = c(1,4))

## 1. Scatter
par(mar = c(4,4,1,1))
plot(cndf[,j] ~ cndf[,i],
     pch = 20, col = coli, las = 1,
     xlab = paste0(i, " CN"), ylab = paste0(j, " CN"),
     xlim = c(0, max_val), ylim = c(0, max_val),
     frame = FALSE)

lower_j <- quantile(cndf[,j], 0.025)
upper_j <- quantile(cndf[,j], 0.975)
lower_i <- quantile(cndf[,i], 0.025)
upper_i <- quantile(cndf[,i], 0.975)
df <- cndf[cndf[,j] >= lower_j & cndf[,j] <= upper_j & cndf[,i] >= lower_i & cndf[,i] <= upper_i,]

fit <- lm(df[,j] ~ df[,i])
abline(fit, lwd = 2/3, lty = 2)
ct <- cor.test(df[,i], df[,j])
text(max_val/20, max_val,
     sprintf("y = %.2f + %.2f x\nr = %.2f, p = %.2g",
             coef(fit)[1], coef(fit)[2],
             ct$estimate, ct$p.value),
     adj = 0)

## 2. Top marginal histogram
par(mar = c(0,4,2,1))
hi <- hist(cndf[,i], breaks = seq(0, max_val, max_val/50), plot = FALSE)
barplot(hi$counts,
        space = 0,
        border = NA,
        axes = FALSE,
        ylim = c(0, max(hi$counts)),
        main = paste0(w, ", ", nrow(cndf), " cells"),
        col = coli)
axis(2, las = 1)

## 3. Right marginal histogram
par(mar = c(4,1,1,1))
hj <- hist(cndf[,j], breaks = seq(0, max_val, max_val/50), plot = FALSE)
barplot(hj$counts,
        horiz = TRUE,
        space = 0,
        border = NA,
        axes = FALSE,
        xlim = c(0, max(hj$counts)),
        main = "",
        col = coli)
axis(1, las = 1)
dev.off()

# Panel B, left
w <- "EL001"

pdf(paste0("Fig3b.left.", w, ".pdf"), height=4, width=4.5)
layout(matrix(1:2, nrow=2), heights=c(0.4, 0.6))

## 1. Top pseudobulk CN plot
qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)
tempdf <- read.csv(paste0("./normalization/subclonal_cn_v4/", w, "_clusters_0.7_1000_subclonal_cn.gc_corrected.csv"))

par(mar = c(0.5, 4.1, 4.1, 2.1))
plot(NULL, xlim = c(hg19$add_values[1], hg19$add_values[nrow(hg19)]+hg19$length[nrow(hg19)]), ylim = c(0, 30), ylab = "Copy number", xlab = "", las=1, xaxt = "n", frame = FALSE, main = paste0(w, "; # hq cancer cells = ", indi$final_cancer_cell_ct[indi$individual == w]))
axis(side = 1, labels = F, at = hg19$add_values)
#text(x = hg19$text_location[1:24], y = -3, labels = hg19$chr[1:24], xpd=NA)
abline(v = hg19$add_values, col = "gray")
tempdf$state[tempdf$state < 0] <- 0
tempdf$state[tempdf$state >= 11] <- 11
segments(x0 = tempdf$plot.start[tempdf$copy != "NaN" & tempdf$copy < 30], y0 = tempdf$copy[tempdf$copy != "NaN" & tempdf$copy < 30], x1 = tempdf$plot.end[tempdf$copy != "NaN" & tempdf$copy < 30], y1 = tempdf$copy[tempdf$copy != "NaN" & tempdf$copy < 30], col = modcol[tempdf$state[tempdf$copy != "NaN" & tempdf$copy < 30]+1], lwd = 1)
if (nrow(tempdf[tempdf$copy != "NaN" & tempdf$copy >= 30,]) != 0){
  segments(x0 = tempdf$plot.start[tempdf$copy != "NaN" & tempdf$copy >= 30], y0 = 30, x1 = tempdf$plot.end[tempdf$copy != "NaN" & tempdf$copy >= 30], y1 = 30, col = "navy", lwd = 1)
}

## 2. Bottom single-cell CN heatmap 
df_final <- fread(paste0("./normalization/cn_dataframe/", w, ".normalized.gc_corrected.filters.annotated.csv"))
df_final$cell_id <- str_sub(df_final$cell_id, -15)

df_final <- df_final[
  normal_reads != 0 &
    segment_length >= 1000 &
    centroacro == "no" &
    blacklist == "no" &
    xfilter == "no" &
    cell_id %in% qc$cell_id[qc$decision == "include"]
]

par(mar = c(3, 4.1, 0.1, 2.1))
cell_number <- max(df_final$sort_rank)
plot(NULL, xlim = c(hg19$add_values[1], hg19$add_values[nrow(hg19)]+hg19$length[nrow(hg19)]), ylim = c(0, cell_number), ylab = "Cancer Cell Copy Number", xlab = "", las=1, xaxt = "n", frame = FALSE)
axis(side = 1, labels = F, at = hg19$add_values)
text(x = hg19$text_location[1:24], y = -max(df_final$sort_rank)*0.07, labels = hg19$chr[1:24], xpd=NA)
abline(v = hg19$add_values, col = "gray")
df_final$state[df_final$state < 0] <- 0
df_final$state[df_final$state >= 11] <- 11
rect(df_final$plot.start, cell_number-df_final$sort_rank, df_final$plot.end, cell_number-df_final$sort_rank+1, border = NA, col = modcol[df_final$state+1])
rm(df_final)
dev.off()

# Panel B, right
w <- "EL001"
i <- "EGFR" # x axis
j <- "CCND1" # y axis
coli <- rgb(39/255, 64/255, 139/255, .3)

pdf(paste0("fig3b.right.", w, ".", i, ".", j, ".pdf"), height=4, width=4)
cndf <- read.csv(paste0("./normalization/oncogene_cn_final/", w, ".sv-seg-cn.oncogene.cn.csv"))
cndf$cell_id <- str_sub(cndf$cell_id, -15)
metrics <- read.csv(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
metrics$cell_id <- str_sub(metrics$cell_id, -15)
cndf <- cndf[cndf$cell_id %in% metrics$cell_id[metrics$decision == "include"],]

max_val <- (floor(max(c(max(cndf[,j], cndf[,i])))/50)+1)*50
layout(matrix(c(2,4,1,3), 2, 2, byrow = TRUE), widths = c(4,1), heights = c(1,4))

## 1. Scatter
par(mar = c(4,4,1,1))
plot(cndf[,j] ~ cndf[,i],
     pch = 20, col = coli, las = 1,
     xlab = paste0(i, " CN"), ylab = paste0(j, " CN"),
     xlim = c(0, max_val), ylim = c(0, max_val),
     frame = FALSE)

lower_j <- quantile(cndf[,j], 0.025)
upper_j <- quantile(cndf[,j], 0.975)
lower_i <- quantile(cndf[,i], 0.025)
upper_i <- quantile(cndf[,i], 0.975)
df <- cndf[cndf[,j] >= lower_j & cndf[,j] <= upper_j & cndf[,i] >= lower_i & cndf[,i] <= upper_i,]

fit <- lm(df[,j] ~ df[,i])
abline(fit, lwd = 2/3, lty = 2)
ct <- cor.test(df[,i], df[,j])
text(max_val/20, max_val,
     sprintf("y = %.2f + %.2f x\nr = %.2f, p = %.2g",
             coef(fit)[1], coef(fit)[2],
             ct$estimate, ct$p.value),
     adj = 0)

## 2. Top marginal histogram
par(mar = c(0,4,2,1))
hi <- hist(cndf[,i], breaks = seq(0, max_val, max_val/50), plot = FALSE)
barplot(hi$counts,
        space = 0,
        border = NA,
        axes = FALSE,
        ylim = c(0, max(hi$counts)),
        main = paste0(w, ", ", nrow(cndf), " cells"),
        col = coli)
axis(2, las = 1)

## 3. Right marginal histogram
par(mar = c(4,1,1,1))
hj <- hist(cndf[,j], breaks = seq(0, max_val, max_val/50), plot = FALSE)
barplot(hj$counts,
        horiz = TRUE,
        space = 0,
        border = NA,
        axes = FALSE,
        xlim = c(0, max(hj$counts)),
        main = "",
        col = coli)
axis(1, las = 1)
dev.off()

# Panel C
w <- "GBM0721"
i <- "128632A-R14-C23"
z <- 4
regions <- "1-201000000-207000000;7-52000000-59000000"
geneofi <- c("MDM4", "EGFR")

svsketch.dlp.sv.seg.cn.onecell(w, i, z, regions, geneofi) # This single-cell plotter function is available in Function part at the end of this script.

# Panel D
w <- "EL001"
i <- "130113A-R55-C12"
z <- 6
regions <- "7-52000000-59000000;11-67000000-72000000;12-55000000-60000000;12-64000000-70000000"
geneofi <- c("EGFR", "CCND1", "CDK4", "MDM2")

svsketch.dlp.sv.seg.cn.onecell(w, i, z, regions, geneofi) # This single-cell plotter function is available in Function part at the end of this script.

# Panel E
w <- "EL001"
qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)

cnvdf <- read.csv(paste0("./normalization/region_cn_final/", w, ".sv-seg-cn.region.cn.csv"))
cnvdf$cell_id <- str_sub(cnvdf$cell_id, -15)
cnvdf <- cnvdf[cnvdf$cell_id %in% qc$cell_id[qc$decision == "include"],]
cnvdf <- cnvdf[,colnames(cnvdf)[colnames(cnvdf) != "chr8_32000001_40000000"]] # excluding non-ecDNA segment

tempdf <- as.data.frame(colnames(cnvdf)[grepl("chr", colnames(cnvdf))])
colnames(tempdf) <- "id"
tempdf$chr_num <- NA
tempdf$start <- NA
tempdf$end <- NA
for (i in 1:nrow(tempdf)){
  tempdf$chr_num[i] <- as.numeric(ifelse(gsub("chr", "", strsplit(tempdf$id[i], "_", fixed = T)[[1]][1]) == "X", "23", gsub("chr", "", strsplit(tempdf$id[i], "_", fixed = T)[[1]][1])))
  tempdf$start[i] <- as.numeric(strsplit(tempdf$id[i], "_", fixed = T)[[1]][2])
  tempdf$end[i] <- as.numeric(strsplit(tempdf$id[i], "_", fixed = T)[[1]][3])
}
tempdf <- tempdf[order(tempdf$start, decreasing = F),]
tempdf <- tempdf[order(tempdf$chr_num, decreasing = F),]
infoline <- tempdf$id

tempdf <- read.csv(paste0("./pairwise_cn/", w, ".paired.allregions_final.csv"))

mat <- matrix(ncol = length(infoline), nrow = length(infoline))
for (k in 1:nrow(tempdf)){
  mat[which(tempdf$region2[k] == infoline),length(infoline)+1-which(tempdf$region1[k] == infoline)] <- ifelse(tempdf$region2[k] == tempdf$region1[k], 1, ifelse(tempdf$pearson_pval[k] < 0.05, tempdf$pearson_r[k], 0))  
}
colnames(mat) <- rev(infoline)
rownames(mat) <- infoline

ht <- Heatmap(
  mat,
  col = colorRamp2(
    breaks = c(0, max(mat, na.rm = T) / 2, max(mat, na.rm = T)),
    colors = viridis::viridis_pal(option = "C")(3) 
  ),
  column_title = "Correlation, clustering",
  name = "Pearson's cor r",
  row_names_side = "left",
  row_names_gp = grid::gpar(fontsize = 6),
  column_names_gp = grid::gpar(fontsize = 6)
)

ht <- draw(ht)
infoline <- column_order(ht)

pdf(paste0("Fig3e.", w, ".pdf"), height=4.2, width=5)
ht
dev.off()

# Panel F, left
w <- "EL001"
geneofi <- c("EGFR", "CCND1", "CDK4", "MDM2")
qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)

celldf <- read.csv("./for-jake/EL001_ordered_loadings.csv")
colnames(celldf) <- c("cell_id", "e1", "e2")
for (i in rev(c("e1", "e2"))){
  celldf <- celldf[order(celldf[,i], decreasing = T),]
}
celldf <- celldf[celldf$cell_id %in% qc$cell_id[qc$decision == "include"],]
celldf$sort_rank <- c(nrow(celldf):1)
celldf <- celldf[order(celldf$sort_rank, decreasing = T),]

df <- east[east$individual == w & east$chromosome %in% c("7", "11", "12"),]
setDT(df)
colnames(df)[2] <- "chr"
setDT(celldf)

df_final <- fread(paste0("./normalization/cn_dataframe/", w, ".normalized.gc_corrected.filters.annotated.csv"))
df_final$cell_id <- str_sub(df_final$cell_id, -15)
df_final <- df_final[cell_id %in% celldf$cell_id]
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

pdf(paste0("Fig3f.left.", w, ".pdf"), width = 5, height = 3.5)
cell_number <- max(df_final$sort_rank)
plot(NULL, xlim = c(0, df$add_value[nrow(df)] + df$seglen[nrow(df)]), ylim = c(0, cell_number), ylab = "Cancer Cell Copy Number", xlab = "", las=1, xaxt = "n", frame = FALSE)
text(x = hg19$text_location[1:24], y = -max(df_final$sort_rank)*0.07, labels = hg19$chr[1:24], xpd=NA)
df_final$state[df_final$state < 0] <- 0
df_final$logcopy <- log10(df_final$state + 1)
df_final$linearcopy <- df_final$copy
df_final$linearcopy[df_final$linearcopy > 500] <- 500

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

pdf(paste0("Fig3f.right.", w, ".pdf"), width = 4, height = 2)
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

# Panel F, right
mat <- as.matrix(celldf[, c(2:3)])
col_names <- colnames(mat) 
min_val <- min(mat)
max_val <- 200

col_list <- list(
  "e1" = circlize::colorRamp2(c(min_val, max_val), c("white", "#2378b5")), # Blue
  "e2" = circlize::colorRamp2(c(min_val, max_val), c("white", "#f57f20"))  # Orange
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

pdf(paste0("Fig3f.right.", w, ".pdf"), width = 2, height = 3.5)
draw(h_list)
dev.off()

# Panel G
w <- "Lx516"
chrs <- c("3", "11", "13", "22")

k <- 23
qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)

df_final <- fread(paste0("./normalization/cn_dataframe/", w, ".normalized.gc_corrected.filters.annotated.csv"))
df_final$cell_id <- str_sub(df_final$cell_id, -15)
df_final <- df_final[cell_id %in% qc$cell_id]

cl <- fread(paste0("./results/sv_gc_normalized_oct_18/", w, "/tables/", w, "_clusters_0.7_1000.csv"))
cl$cell_id <- str_sub(cl$cell_id, -15)
cl <- cl[cl$cell_id %in% qc$cell_id[qc$decision == "include"] & cl$cluster_label != "F",]

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

isv <- read.csv(paste0("./hmf/gridss_somatic/", w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.annotated.bedpe"), header=T, as.is=T, sep="\t")
isv <- isv[!(isv$type == "TRA" & isv$vf <= 4),]
isv <- isv[!(isv$type %in% c("t2tFBI", "h2hFBI") & isv$vf <= 4),]
colnames(isv)[colnames(isv) == "ori1"] <- "strand1"
colnames(isv)[colnames(isv) == "ori2"] <- "strand2"
isv$pos1 <- (isv$start1 + isv$end1)/2
isv$pos2 <- (isv$start2 + isv$end2)/2

isv <- isv[isv$chr1 %in% chrs | isv$chr2 %in% chrs,]
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

pdf(paste0("Fig3g.", w, "_", paste(chrs, collapse = "-"), ".pdf"), width=8.5, height=4)
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

par(new = T)
if (k < 5){
  k <- 5
}
plot(0, type='n', xlim=c(0, hg19$length[hg19$chr == deltadf$chrs[nrow(deltadf)]] + deltadf$delta[nrow(deltadf)]), ylim=c(0,1.1*k), frame=F, xaxt="n", las=1, xlab=paste0("Position on chromosomes ", paste(chrs, collapse = ", ")), ylab="Total copy number")

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

df <- df_final[cell_id %in% cl$cell_id & chr %in% chrs]
cnvdf <- consensuscopynumber(df)
colnames(cnvdf) <- c("chromosome", "start", "end", "cell_id", "state", "totalCopyNumber")

for (j in 1:nrow(cnvdf)){
  cnvdf$start[j] <- cnvdf$start[j] + deltadf$delta[deltadf$chrs == cnvdf$chromosome[j]]
  cnvdf$end[j] <- cnvdf$end[j] + deltadf$delta[deltadf$chrs == cnvdf$chromosome[j]]
}
cnvdf <- cnvdf[order(cnvdf$start, decreasing = F),]
cnvdf$seglen <- cnvdf$end - cnvdf$start
cnvdf <- cnvdf[cnvdf$seglen >= 500,]
cnvdf$state[cnvdf$state > 10] <- 11
segments(cnvdf$start[cnvdf$start != cnvdf$end], cnvdf$totalCopyNumber[cnvdf$start != cnvdf$end], cnvdf$end[cnvdf$start != cnvdf$end], cnvdf$totalCopyNumber[cnvdf$start != cnvdf$end], col=modcol[cnvdf$state + 1], lwd=1, lend=0)

for (j in c("CCND1", "KDM2A", "KMT2A")){
  if (genedf$chromosome[genedf$gene == j] %in% chrs){
    segments(x0 = genedf$start[genedf$gene == j] + deltadf$delta[deltadf$chrs == genedf$chromosome[genedf$gene == j]], y0 = 1.1*k, x1 = genedf$end[genedf$gene == j] + deltadf$delta[deltadf$chrs == genedf$chromosome[genedf$gene == j]], y1 = 1.1*k, col = "darkgreen", lwd = 2, xpd = T)
    text(x = genedf$start[genedf$gene == j] + deltadf$delta[deltadf$chrs == genedf$chromosome[genedf$gene == j]], y = 1.2*k, labels = j, xpd = T)
  }
}

for (j in chrs){
  rect(cytoband$start[cytoband$chr == j] + deltadf$delta[deltadf$chrs == j], (-0.22)*k, cytoband$end[cytoband$chr == j] + deltadf$delta[deltadf$chrs == j], (-0.17)*k, col = cytoband$color[cytoband$chr == j], xpd = T, lwd = 0.5)
}
rm(df_final)
rm(df)
dev.off()


# Panel H
w <- "Lx516"

cl <- fread(paste0("./results/sv_gc_normalized_oct_18/", w, "/tables/", w, "_clusters_0.7_1000.csv"))
cl$cell_id <- str_sub(cl$cell_id, -15)
qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)
cl <- cl[cell_id %in% qc$cell_id[qc$decision == "include"]]
cnvdf <- read.csv(paste0("./normalization/region_cn_final/", w, ".sv-seg-cn.region.cn.csv"))
cnvdf$cell_id <- str_sub(cnvdf$cell_id, -15)
cnvdf <- cnvdf[!(cnvdf$cell_id %in% cl$cell_id[cl$cluster_label == "F"]),] # Excluding clonally unrelated cells
cnvdf <- cnvdf[cnvdf$cell_id %in% qc$cell_id[qc$decision == "include"],]

tempdf <- as.data.frame(colnames(cnvdf)[grepl("chr", colnames(cnvdf))])
colnames(tempdf) <- "id"
tempdf$chr_num <- NA
tempdf$start <- NA
tempdf$end <- NA
for (i in 1:nrow(tempdf)){
  tempdf$chr_num[i] <- as.numeric(ifelse(gsub("chr", "", strsplit(tempdf$id[i], "_", fixed = T)[[1]][1]) == "X", "23", gsub("chr", "", strsplit(tempdf$id[i], "_", fixed = T)[[1]][1])))
  tempdf$start[i] <- as.numeric(strsplit(tempdf$id[i], "_", fixed = T)[[1]][2])
  tempdf$end[i] <- as.numeric(strsplit(tempdf$id[i], "_", fixed = T)[[1]][3])
}
tempdf <- tempdf[order(tempdf$start, decreasing = F),]
tempdf <- tempdf[order(tempdf$chr_num, decreasing = F),]
infoline <- tempdf$id

tempdf <- read.csv(paste0("./pairwise_cn/", w, ".paired.allregions_final.csv"))
mat <- matrix(ncol = length(infoline), nrow = length(infoline))
for (k in 1:nrow(tempdf)){
  mat[which(tempdf$region2[k] == infoline),length(infoline)+1-which(tempdf$region1[k] == infoline)] <- ifelse(tempdf$region2[k] == tempdf$region1[k], 1, ifelse(tempdf$pearson_pval[k] < 0.05, tempdf$pearson_r[k], 0))  
}
colnames(mat) <- rev(infoline)
rownames(mat) <- infoline

ht <- Heatmap(
  mat,
  col = colorRamp2(
    breaks = c(0, max(mat, na.rm = T) / 2, max(mat, na.rm = T)),
    colors = viridis::viridis_pal(option = "C")(3) 
  ),
  column_title = "Correlation, clustering",
  name = "Pearson's cor r", # Changed name for relevance
  row_names_side = "left",
  row_names_gp = grid::gpar(fontsize = 6),
  column_names_gp = grid::gpar(fontsize = 6)
)

ht <- draw(ht)
infoline <- column_order(ht)

pdf(paste0("Fig3h.", w, ".pdf"), height=4.2, width=5)
ht
dev.off()

# Panel I: conceptual illustration

# Function: single-cell plotter, for one cell
svsketch.dlp.sv.seg.cn.onecell <- function(w, i, z, regions, geneofi){
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
  
  qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
  qc$cell_id <- str_sub(qc$cell_id, -15)
  
  df_final <- fread(paste0("./normalization/cn_dataframe/", w, ".normalized.gc_corrected.filters.annotated.csv"))
  df_final$cell_id <- str_sub(df_final$cell_id, -15)
  
  df_final <- df_final[
    normal_reads != 0 &
      segment_length >= 1000 &
      #centroacro == "no" &
      blacklist == "no" &
      xfilter == "no" &
      cell_id %in% qc$cell_id[qc$decision == "include"]
  ]
  
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
  
  for (j in 1:length(regs)){
    df_final$plot.start[df_final$regs == regs[j]] <- df_final$start[df_final$regs == regs[j]] - df$start[df$regs == regs[j]] + df$add_value[df$regs == regs[j]]
    df_final$plot.end[df_final$regs == regs[j]] <- df_final$end[df_final$regs == regs[j]] - df$start[df$regs == regs[j]] + df$add_value[df$regs == regs[j]]
  }
  
  svdf <- fread(paste0("./hmf/gridss_", ifelse(w %in% indi$individual[indi$normal_wgs == "yes"], "somatic/", "tumoronly/"), w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.annotated.bedpe"))
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
    v[is.na(v) & toupper(x) %in% c("M", "MT")] <- 25 
    v[is.na(v)] <- 9999 
    v
  }]
  setorder(svdf, chr_ord, start1, end1)
  svdf[, chr_ord := NULL]
  
  svdf[, reg1 := df[
    .SD, 
    on = .(chr = chr1, start <= pos1, end >= pos1), 
    j = regs, 
    nomatch = NA_character_
  ]]
  
  svdf[, reg2 := df[
    .SD, 
    on = .(chr = chr2, start <= pos2, end >= pos2), 
    j = regs, 
    nomatch = NA_character_
  ]]
  
  filelist <- list.files(
    path = file.path("./sc_sv_genotyping", w),
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
  
  colnames(svmat) <- gsub("_right", "", gsub("_left", "", colnames(svmat)))
  svmat <- svmat[cell_id %in% qc$cell_id[qc$decision == "include"]]
  svmat <- svmat[, .SD, .SDcols = c("cell_id", names(svmat)[names(svmat) %in% unique(svdf$svname)])]
  mat <- as.matrix(svmat[, .SD, .SDcols = 2:ncol(svmat)])
  rownames(mat) <- svmat[["cell_id"]]
  
  filelist <- list.files(
    path = file.path("./sc_sv_genotyping", w),
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
  
  colnames(svmat) <- gsub("_right", "", gsub("_left", "", colnames(svmat)))
  svmat <- svmat[cell_id %in% qc$cell_id[qc$decision == "include"]]
  svmat <- svmat[, .SD, .SDcols = c("cell_id", names(svmat)[names(svmat) %in% unique(svdf$svname)])]
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
  svdf <- svdf[svname %in% colnames(mat)]
  celllist <- qc$cell_id[qc$decision == "include" & qc$cell_id %in% rownames(mat)]
  
  pdf(paste0("Fig3cd.", w, ".", i, ".", paste(chrs, collapse = "-"), ".pdf"), width=z, height=4) # height 4 for Fig. 3 but 3 for Fig. 3c, 3.7 for 4D
  cnvdf <- as.data.frame(df_final[cell_id == i])
  cnvdf <- cnvdf[,c("chr", "plot.start", "plot.end", "copy", "state")]
  k <- max(cnvdf$copy)
  colnames(cnvdf) <- c("chromosome", "plotstart", "plotend", "CopyNumber", "state")
  cnvdf$state[cnvdf$state >= 11] <- 11
  cnvdf$state[cnvdf$state < 0] <- 0
  cnvdf$plotstart <- as.numeric(cnvdf$plotstart)
  cnvdf$plotend <- as.numeric(cnvdf$plotend)
  
  ## Plotting Copy Numbers
  if (k < 5){
    k <- 5
  }
  plot(0, type='n', xlim=c(0, df$add_value[nrow(df)] + df$seglen[nrow(df)]), ylim=c(0,1.1*k), xaxt = 'n', frame=F, las=1, xlab=paste0("Position on chromosome ", paste(unique(df$chr), collapse = "-"), " in Mbp"), ylab="Total copy number", main = paste0(i, " quality = ", qc$quality[qc$cell_id == i]))
  for (j in 1:nrow(df)){
    axis(1, at = c(df$add_value[j], df$add_value[j] + df$seglen[j]), tick = c(df$add_value[j]+1, df$add_value[j] + df$seglen[j]), labels = c(df$start[j], df$end[j]), las=2)
  }
  segments(cnvdf$plotstart, cnvdf$CopyNumber, cnvdf$plotend, cnvdf$CopyNumber, col=modcol[as.numeric(cnvdf$state)+1], lwd=2/3)
  
  ## Genes of Interest
  for (j in geneofi){
    segments(x0 = genedf$start[genedf$gene == j] - df$start[grepl(j, df$gene)] + df$add_value[grepl(j, df$gene)], y0 = 1.1*k, x1 = genedf$end[genedf$gene == j] - df$start[grepl(j, df$gene)] + df$add_value[grepl(j, df$gene)], y1 = 1.1*k, col = "darkblue", lwd = 2, xpd = T)
    text(x = genedf$start[genedf$gene == j] - df$start[grepl(j, df$gene)] + df$add_value[grepl(j, df$gene)], y = 1.2*k, labels = j, xpd = T)
  }
  
  l = max(mat[i,], na.rm = T)
  if (l < 5){
    l <- 5
  }
  
  if (length(colnames(mat)[mat[i,] != 0]) != 0){
    isv <- as.data.frame(svdf[svname %in% colnames(mat)[mat[i,] != 0]])
    isv$pos1 <- as.numeric(isv$pos1)
    isv$pos2 <- as.numeric(isv$pos2)
    isv$plotpos1 <- isv$pos1
    isv$plotpos2 <- isv$pos2
    
    for (j in 1:nrow(df)){
      isv$plotpos1[isv$reg1 == df$regs[j] & !is.na(isv$reg1)] <- isv$plotpos1[isv$reg1 == df$regs[j] & !is.na(isv$reg1)] - df$start[j] + df$add_value[j]
      isv$plotpos2[isv$reg2 == df$regs[j] & !is.na(isv$reg2)] <- isv$plotpos2[isv$reg2 == df$regs[j] & !is.na(isv$reg2)] - df$start[j] + df$add_value[j]
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
          svcol <- rgb(0/255,0/255,255/255,.3)
          x <- isvin$rad[j]*cos(theta)+((isvin$plotpos1[j]+isvin$plotpos2[j])/2)
          y <- l*0.1*sin(theta)+mat[i,isvin$svname[j]]
        } else if (isvin$strand1[j] == "-" & isvin$strand2[j] == "+"){
          svcol <- rgb(0/255,128/255,128/255,.3)
          x <- isvin$rad[j]*cos(theta)+((isvin$plotpos1[j]+isvin$plotpos2[j])/2)
          y <- -l*0.1*sin(theta)+mat[i,isvin$svname[j]]
        } else if (isvin$strand1[j] == "-" & isvin$strand2[j] == "-"){
          svcol <- rgb(220/255,20/255,60/255,.3)
          x <- isvin$rad[j]*cos(theta)+((isvin$plotpos1[j]+isvin$plotpos2[j])/2)
          y <- l*0.1*sin(theta)+mat[i,isvin$svname[j]]
        } else if (isvin$strand1[j] == "+" & isvin$strand2[j] == "+"){
          svcol <- rgb(128/255,128/255,0/255,.3)
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
        svcol <- rgb(85/255,26/255,139/255,.15)
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

