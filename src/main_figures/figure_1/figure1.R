# Figure 1
library(stringr)
library(plotrix)

# Panel A: number of samples by cancer type
indi <- read.csv("dlp.cohort.summary.individual_111925.tsv", header = T, as.is = T, sep = "\t")
table(indi$cancer_type)

# Panel B
df <- as.data.frame(sort(table(indi$cancer_type), decreasing = T))
colnames(df) <- c("cancer_type", "n_patient")
df$cancer_cell_ct <- 0
for (i in 1:nrow(df)) {
  df$cancer_cell_ct[i] <- sum(indi$final_cancer_cell_ct[indi$cancer_type == df$cancer_type[i]])
}
df <- df[order(df$cancer_cell_ct, decreasing = T), ]

pdf("Fig1b.pdf", height = 3.5, width = 5)
par(mar = c(5.1, 8.1, 4.1, 2.1))
barplot(rev(df$cancer_cell_ct), log = "x", las = 1, names.arg = rev(paste0(df$cancer_type, " (n=", df$n_patient, ")")), border = F, xlim = c(1, 100000), col = rev(c("#008B8B", "#CD3278", "#00BFFF", "#BABABA", "#FFB5C5", "#FF9912", "#27408B")), xlab = "Number of high-quality cancer cells", horiz = T, space = 1)
dev.off()

# Panel C
df <- read.csv("modeling/ecdna_model_summary_table_111925.csv")
df <- df[df$gene %in% names(sort(table(df$gene)[table(df$gene) >= 3], decreasing = T))[!(names(sort(table(df$gene)[table(df$gene) >= 3], decreasing = T)) %in% c("ANK1", "NCOA2", "MECOM", "PABPC1", "IKBKB", "JAK3", "CALR", "FGF3", "IDH2"))][1:15], ]
df$gene <- factor(df$gene, levels = names(sort(table(df$gene)[table(df$gene) >= 3], decreasing = T))[!(names(sort(table(df$gene)[table(df$gene) >= 3], decreasing = T)) %in% c("ANK1", "NCOA2", "MECOM", "PABPC1", "IKBKB", "JAK3", "CALR", "FGF3", "IDH2"))][1:15])
df$cancer_type <- factor(df$cancer_type, levels = names(sort(table(df$cancer_type), decreasing = T)))
pdf("Fig1c.pdf", height = 3.5, width = 5)
barplot(t(as.matrix(table(df$gene, df$cancer_type))), las = 2, col = c("#008B8B", "#CD3278", "#27408B", "#BABABA", "#00BFFF", "#FFB5C5", "#FF9912"), border = F, ylab = "Number of patients with HLAMP", space = 1)
dev.off()

# Panel D: conceptual illustration

# Panel E, left
w <- "P-0009535"
j <- "EGFR"

# Panel E, right
w <- "GBM0721"
j <- "EGFR"

# Panel E, common downstream
cndf <- read.csv(paste0("oncogene_cn_final/", w, ".sv-seg-cn.oncogene.cn.csv"))
infoline <- cndf[, j]
infoline <- infoline[!is.na(infoline) & infoline > 0]

pdf(paste0("Fig1d_", w, ".", j, ".pdf"), height = 4, width = 3.5)
hist(infoline, las = 1, breaks = seq(0, ceiling(max(infoline) / 50) * 50, ceiling(max(infoline) / 50)), xlab = "Per cell copy number", main = paste0(w, " (n=", length(infoline), ") ", j), col = "#FF82AB", border = F, prob = T)
lines(density(infoline), lwd = 2, col = "#B0171F")
dev.off()

# Panel F: conceptual illustration

# Panel G
east <- read.csv("modeling/amplified_regions_summary_table_011026.csv")
pdf("Fig1g.pdf", height = 5, width = 5.5)
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), frame = F, las = 1, xlab = "Mass in window", ylab = "ecDNA score")
for (i in which(east$mean_cn > 6)) {
  draw.circle(x = east$mass_in_window[i], y = east$score[i], radius = log(east$n_cells[i]) / 500, col = c("#008B8B77", "#CD327877", "#00BFFF77", "#FFB5C577", "#BABABA77", "#27408B77", "#FF991277")[east$cancer_type[i]], border = F)
}
abline(b = 0.859, a = -0.184, lwd = 2 / 3) # Decision boundary line: y = 0.859x + -0.184

## Legend
draw.circle(x = 1, y = 1, radius = log(30) / 500, col = "#BABABAAA", border = F)
draw.circle(x = 1, y = 0.95, radius = log(100) / 500, col = "#BABABAAA", border = F)
draw.circle(x = 1, y = 0.9, radius = log(300) / 500, col = "#BABABAAA", border = F)
draw.circle(x = 1, y = 0.85, radius = log(1000) / 500, col = "#BABABAAA", border = F)
draw.circle(x = 1, y = 0.8, radius = log(3000) / 500, col = "#BABABAAA", border = F)
draw.circle(x = 1, y = 0.75, radius = log(10000) / 500, col = "#BABABAAA", border = F)
dev.off()

# Panel H
df <- east[east$mean_cn > 6, ]
tempdf <- as.data.frame(sort(table(df$individual), decreasing = T))
colnames(tempdf) <- c("individual", "n_amp")
tempdf$n_07subclonal_ICamp <- 0
tempdf$n_07clonal_ICamp <- 0
tempdf$n_07subclonal_ecDNA <- 0
tempdf$n_07clonal_ecDNA <- 0
tempdf$cancer_type <- ""

for (i in 1:nrow(tempdf)) {
  print(paste0("Processing ", tempdf$individual[i]))
  tempdf$cancer_type[i] <- indi$cancer_type[indi$individual == tempdf$individual[i]]
  cl <- read.csv(paste0("normalization/acgr_subclonal_cn_v4/", tempdf$individual[i], "_clusters_0.7_1000_subclonal_cn_acgr.csv"))
  clonelist <- colnames(cl)[grepl("clone_", colnames(cl))]
  df <- east[east$individual == tempdf$individual[i] & east$mean_cn >= 6, ]
  if (tempdf$individual[i] != "2765_2") {
    cndf <- read.csv(paste0("sv_gc_normalized_oct_18/", tempdf$individual[i], "/tables/", tempdf$individual[i], "_clusters_0.7_1000.csv"))
  } else {
    cndf <- read.csv("sv_cna_clones_nov_14_v2/2765_2/tables/2765_2_clusters_0.7_1000.csv")
  }
  cndf$cell_id <- str_sub(cndf$cell_id, -15)
  qc <- read.csv(paste0("metadata/doublet_invariant_filtered/", tempdf$individual[i], "_qc.filtered.csv"))
  qc$cell_id <- str_sub(qc$cell_id, -15)
  cndf <- cndf[cndf$cell_id %in% qc$cell_id[qc$decision == "include"], ]
  for (j in 1:nrow(df)) {
    chri <- df$chromosome[j]
    starti <- df$start[j]
    endi <- df$end[j]
    infoline <- sapply(cl[cl$chr == chri & cl$start >= starti & cl$end <= endi, clonelist], function(x) {
      weighted.mean(x, w <- cl$seglen[cl$chr == chri & cl$start >= starti & cl$end <= endi], na.rm = T)
    })
    if (df$decision[j] == "ICamp" & nrow(cndf[cndf$cluster_label %in% gsub("clone_", "", names(infoline)[infoline > 6]), ]) / nrow(cndf[cndf$cluster_label %in% gsub("clone_", "", clonelist), ]) >= 0.9) {
      tempdf$n_07clonal_ICamp[i] <- tempdf$n_07clonal_ICamp[i] + 1
    } else if (df$decision[j] == "ICamp" & nrow(cndf[cndf$cluster_label %in% gsub("clone_", "", names(infoline)[infoline > 6]), ]) / nrow(cndf[cndf$cluster_label %in% gsub("clone_", "", clonelist), ]) < 0.9) {
      tempdf$n_07subclonal_ICamp[i] <- tempdf$n_07subclonal_ICamp[i] + 1
    } else if (df$decision[j] == "ecDNA" & nrow(cndf[cndf$cluster_label %in% gsub("clone_", "", names(infoline)[infoline > 6]), ]) / nrow(cndf[cndf$cluster_label %in% gsub("clone_", "", clonelist), ]) >= 0.9) {
      tempdf$n_07clonal_ecDNA[i] <- tempdf$n_07clonal_ecDNA[i] + 1
    } else if (df$decision[j] == "ecDNA" & nrow(cndf[cndf$cluster_label %in% gsub("clone_", "", names(infoline)[infoline > 6]), ]) / nrow(cndf[cndf$cluster_label %in% gsub("clone_", "", clonelist), ]) < 0.9) {
      tempdf$n_07subclonal_ecDNA[i] <- tempdf$n_07subclonal_ecDNA[i] + 1
    }
    rm("chri", "starti", "endi", "infoline")
  }
}
mat <- as.matrix(tempdf[, c(3:6)])
row.names(mat) <- tempdf$individual
tempdf$cancer_type <- factor(tempdf$cancer_type, levels = c("HGSOC", "HR-BC", "EGFR+LUAD", "EXP", "HR+BC", "SCLC", "GBM"))

pdf("Fig1h.upper.pdf", height = 4, width = 6)
barplot(t(mat), las = 2, col = c("#FFBBFF", "#8B008B", "#F0E68C", "#FF9912"), border = NA, ylab = "Number of amplified regions")
legend("topright", rev(colnames(mat)), fill = rev(c("#FFBBFF", "#8B008B", "#F0E68C", "#FF9912")))
dev.off()

pdf("Fig1h.lower.pdf", height = 3, width = 6)
image(as.matrix(as.numeric(tempdf$cancer_type)), col = c("#008B8B", "#CD3278", "#27408B", "#BABABA", "#FFB5C5", "#00BFFF", "#FF9912"), axes = F)
legend(
  x = "topright",
  inset = c(-0.15, 0),
  legend = paste(names(table(tempdf$cancer_type))),
  fill = c("#008B8B", "#CD3278", "#27408B", "#BABABA", "#FFB5C5", "#00BFFF", "#FF9912"),
  title = "Cancer type",
  bg = "transparent",
  xpd = TRUE
)
dev.off()
