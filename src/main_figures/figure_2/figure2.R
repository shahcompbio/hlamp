# Figure 2
library(stringr)
indi <- read.csv("./metadata/dlp.cohort.summary.individual_111925.tsv", header=T, as.is=T, sep="\t")
east <- read.csv("./modeling/amplified_regions_summary_table_011026.csv")

# Panel A
pdf("Fig2a.pdf", height=4, width=3.5)
bp <- barplot(table(east$n_peaks[east$mean_cn > 6 & east$decision == "ICamp" & east$n_peaks != 0]), las = 1, col = "#8B1C62", border = NA, main = "Number of peaks in ICamps (mean CN > 6)", space = 1)
axis(4, at = c(0, 0.1, 0.2, 0.3) * sum(table(east$n_peaks[east$mean_cn > 6 & east$decision == "ICamp" & east$n_peaks != 0])), labels = c("0.0", "0.1", "0.2", "0.3"), las = 1)
mtext("Fraction", side = 4, line = 3)
dev.off()

# Panel B
w <- "SPECTRUM-OV-044"
j <- "ERBB2"

cndf <- read.csv(paste0("./normalization/oncogene_cn_final/", w, ".sv-seg-cn.oncogene.cn.csv"))
infoline <- cndf[,j]
infoline <- infoline[!is.na(infoline) & infoline >0]

pdf(paste0("Fig2b_", w, ".", j, ".pdf"), height=4, width=3.5)
hist(infoline, las=1, breaks = seq(0,ceiling(max(infoline)/50)*50,ceiling(max(infoline)/50)), xlab = "Per cell copy number", main = paste0(w, " (n=", length(infoline), ") ", j), col = "#FF82AB", border = F, prob = T)
lines(density(infoline), lwd=2, col="#B0171F")
dev.off()

# Panel C
df <- indi
df$sd_TP53 <- NA
df$sd_RB1 <- NA
for (j in 1:nrow(df)){
  w <- df$individual[j]
  print(paste0("Processing ", w))
  cndf <- read.csv(paste0("./normalization/ref_gene_cn_v4/", w, ".sv-seg-cn.ref_gene.cn.csv"))
  for (i in c("TP53", "RB1")){
    df[j,paste0("sd_", i)] <- sd(cndf[,i])
  }
}

options(scipen = 999)
tempdf <- east[east$mean_cn > 6,c("sd_cn", "process")]
tempdf$process <- as.character(tempdf$process)
infodf <- as.data.frame(c(df$sd_TP53, df$sd_RB1))
colnames(infodf) <- "sd_cn"
infodf$process <- c(rep("TP53", nrow(df)), rep("RB1", nrow(df)))
tempdf <- rbind(tempdf, infodf)
tempdf$process <- factor(tempdf$process, levels = rev(c("ecDNA", "CTX", "BFB", "CGR-NOS", "CGR-TRA", "TD", "UNK", "TP53", "RB1")))

pdf("Fig2c.pdf", height=4, width=5)
par(mar = c(5.1, 11.1, 2.1, 2.1))
boxplot(tempdf$sd_cn ~ tempdf$process, log = "x", frame = F, las = 1, outline = F, ylab = "", xlab = "Standard deviation of copy number", horizontal = T, ylim = c(0.1, 150), names = rev(c(paste0("ecDNA (", nrow(tempdf[tempdf$process == "ecDNA",]), ")"), paste0("chromothripsis (", nrow(tempdf[tempdf$process == "CTX",]), ")"), paste0("BFB (", nrow(tempdf[tempdf$process == "BFB",]), ")"), paste0("CGR-NOS (", nrow(tempdf[tempdf$process == "CGR-NOS",]), ")"), paste0("CGR-TRA (", nrow(tempdf[tempdf$process == "CGR-TRA",]), ")"), paste0("tandem duplication (", nrow(tempdf[tempdf$process == "TD",]), ")"), paste0("unknown (", nrow(tempdf[tempdf$process == "UNK",]), ")"), "TP53 (112)", "RB1 (112)")), col = c("#C0C0C0", "#808080", "#8B4C39", "#CD7054", "#EE8262", "#FF8C69", "#DC143C", "#800080", "#FF9912"), lwd = 1/3)
stripchart(tempdf$sd_cn ~ tempdf$process, pch = 20, col = rgb(0,0,0,.1), method = "jitter", vertical = F, add = T, log = "x")
dev.off()

# Panel D: conceptual illustration

# Panel E
df <- read.csv("./tables/skyline_jake_x_clustered.csv")
df <- df[!is.na(df$x_cluster),]
df$cluster <- "numeric"
df$cluster[df$x_cluster == 1] <- "structural"
df$cluster[df$x_cluster_prob <= 0.95] <- "ambiguous"
table(df$cluster)
table(df$pair[df$cluster == "numeric"])

pdf("Fig2e", height=3.8, width=5)
plot(df$x, df$y, ylim = c(0,250), xlab = "Delta cosine similarity for HLAMP region", ylab = "Copy number difference at HLAMP peak (at least 2)", las = 1, frame = F, main = "HLAMP CN comparison between clones", pch = 20, col = c(rgb(209/255, 224/255, 122/255, 0.3), rgb(200/255, 200/255, 200/255, 0.3), rgb(247/255, 175/255, 156/255, 0.3))[factor(df$cluster, levels = c("numeric", "ambiguous", "structural"))])
legend("topleft", legend = c(paste0("Numeric (n=", nrow(df[df$cluster == "numeric",]), ")"), paste0("Intermediate (n=", nrow(df[df$cluster == "ambiguous",]), ")"), paste0("Structural (n=", nrow(df[df$cluster == "structural",]), ")")), pch = 20, col = c(rgb(209/255, 224/255, 122/255, 0.3), rgb(200/255, 200/255, 200/255, 0.3), rgb(247/255, 175/255, 156/255, 0.3)), box.lty = 0)
dev.off()

# Panel F
w <- "EL011"
clonelist <- c("C", "H")
chrs <- c("19")

# Panel G
w <- "SPECTRUM-OV-002"
clonelist <- c("A", "B")
chrs <- c("17")

# Common code for Panel F/G generation
pdf(paste0("fig2fg.", w, ".pdf"), height=3.8, width=5)
qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)

df_final <- fread(paste0("./normalization/cn_dataframe/", w, ".normalized.gc_corrected.filters.annotated.csv"))
df_final <- df_final[
  normal_reads != 0 &
    segment_length >= 1000 &
    centroacro == "no" &
    blacklist == "no" &
    xfilter == "no"
]
df_final$cell_id <- str_sub(df_final$cell_id, -15)
df_final <- df_final[cell_id %in% qc$cell_id]

cl <- read.csv(paste0("./results/sv_gc_normalized_oct_18/", w, "/tables/", w, "_clusters_0.7_1000.csv"))
cl$cell_id <- str_sub(cl$cell_id, -15)
colnames(cl)[2] <- "clone_id"
cl$clone_id[!(cl$clone_id %in% clonelist)] <- "0"

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

isv <- read.csv(paste0("./hmf/gridss_", ifelse(indi$normal_wgs[indi$individual == w] == "yes", "somatic", "tumoronly"), "/", w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.annotated.bedpe"), header=T, as.is=T, sep="\t")
isv <- isv[!(isv$type == "TRA" & isv$vf <= 4),]
isv <- isv[!(isv$type %in% c("t2tFBI", "h2hFBI") & isv$vf <= 4),]
colnames(isv)[colnames(isv) == "ori1"] <- "strand1"
colnames(isv)[colnames(isv) == "ori2"] <- "strand2"
isv$pos1 <- (isv$start1 + isv$end1)/2
isv$pos2 <- (isv$start2 + isv$end2)/2

isv <- isv[isv$chr1 %in% chrs | isv$chr2 %in% chrs,]
if (nrow(isv) != 0){
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
} else {
  max_val <- 5
}

plot(0, type='n', xlim=c(0, hg19$length[hg19$chr == deltadf$chrs[nrow(deltadf)]] + deltadf$delta[nrow(deltadf)]), ylim=c(0,1.1*max_val), frame=F, axes = F, xlab = "", ylab = "")
axis(4, las=1)

if (nrow(isv) != 0){
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
  
  if (nrow(isvin) != 0){
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
  }
}

par(new = T)

k <- 5
for (l in clonelist){
  df <- df_final[cell_id %in% cl$cell_id[cl$clone_id == l] & chr %in% chrs]
  cnvdf <- consensuscopynumber(df)
  if (k < max(cnvdf$copy, na.rm = T)){
    k <- max(cnvdf$copy, na.rm = T)
  }
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

df <- df_final[cell_id %in% cl$cell_id[cl$clone_id == clonelist[1]] & chr %in% chrs]
cnvdf <- consensuscopynumber(df)
colnames(cnvdf) <- c("chromosome", "start", "end", "cell_id", "state", "copy1")
cnvdf$copy2 <- consensuscopynumber(df_final[cell_id %in% cl$cell_id[cl$clone_id == clonelist[2]] & chr %in% chrs])$copy

for (j in 1:nrow(cnvdf)){
  cnvdf$start[j] <- cnvdf$start[j] + deltadf$delta[deltadf$chrs == cnvdf$chromosome[j]]
  cnvdf$end[j] <- cnvdf$end[j] + deltadf$delta[deltadf$chrs == cnvdf$chromosome[j]]
}

cnvdf$seglen <- cnvdf$end - cnvdf$start
cnvdf <- cnvdf[order(cnvdf$start, decreasing = F),]
cnvdf <- cnvdf[cnvdf$seglen >= 500,]

segments(cnvdf$start, cnvdf$copy2, cnvdf$end, cnvdf$copy2, col=base_colors[1], lwd=1, lend=0)
for (j in 2:nrow(cnvdf)){
  segments(cnvdf$start[j], cnvdf$copy2[j-1], cnvdf$start[j], cnvdf$copy2[j], col=base_colors[1], lwd=1, lend=0)
}
segments(cnvdf$start, cnvdf$copy1, cnvdf$end, cnvdf$copy1, col=base_colors[3], lwd=1, lend=0)
for (j in 2:nrow(cnvdf)){
  segments(cnvdf$start[j], cnvdf$copy1[j-1], cnvdf$start[j], cnvdf$copy1[j], col=base_colors[3], lwd=1, lend=0)
}
rect(cnvdf$start, cnvdf$copy1, cnvdf$end, cnvdf$copy2, col = "#ab332922", border = NA)
rect(cnvdf$start[cnvdf$copy1 > 1], 1, cnvdf$end[cnvdf$copy1 > 1], cnvdf$copy1[cnvdf$copy1 > 1], col = "#00688B22", border = NA)

legend("topleft", legend = paste0("clone ", clonelist), col = base_colors[c(1,3)], lty = 1, bty = "n", border = NA)

# GENES OF INTEREST
for (j in c("CCNE1", "ERBB2")){
  if (genedf$chromosome[genedf$gene == j] %in% chrs){
    segments(x0 = genedf$start[genedf$gene == j] + deltadf$delta[deltadf$chrs == genedf$chromosome[genedf$gene == j]], y0 = 1.1*k, x1 = genedf$end[genedf$gene == j] + deltadf$delta[deltadf$chrs == genedf$chromosome[genedf$gene == j]], y1 = 1.1*k, col = "darkgreen", lwd = 2, xpd = T)
    text(x = genedf$start[genedf$gene == j] + deltadf$delta[deltadf$chrs == genedf$chromosome[genedf$gene == j]], y = 1.2*k, labels = j, xpd = T)
  }
}

# IDEOGRAM
for (j in chrs){
  rect(cytoband$start[cytoband$chr == j] + deltadf$delta[deltadf$chrs == j], (-0.22)*k, cytoband$end[cytoband$chr == j] + deltadf$delta[deltadf$chrs == j], (-0.17)*k, col = cytoband$color[cytoband$chr == j], xpd = T, lwd = 0.5)
}
rm(df_final)
rm(df)
dev.off()

# Panel H, left -- we renamed subclones in the manuscript to follow the flow.
w <- "SPECTRUM-OV-052"
j <- "CCND1"
clonelist <- c("C", "A", "B", "0")

## Loading the data
cndf <- read.csv(paste0("./normalization/oncogene_cn_final/", w, ".sv-seg-cn.oncogene.cn.csv"))
cndf$cell_id <- str_sub(cndf$cell_id, -15)
metrics <- read.csv(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
metrics$cell_id <- str_sub(metrics$cell_id, -15)
cndf <- cndf[cndf$cell_id %in% metrics$cell_id[metrics$decision == "include"],]

cl <- read.csv(paste0("./tables/", w, "/", w, "_clusters_0.3_1000.csv"))
cl$cell_id <- str_sub(cl$cell_id, -15)
colnames(cl)[2] <- "clone_id"
cl$clone_id[cl$clone_id %in% c("C", "D", "E", "F", "H", "I")] <- "0"
cl$clone_id[cl$clone_id == "G"] <- "C"
cndf$clone_id <- ""
for (i in which(cndf$cell_id %in% cl$cell_id)){
  cndf$clone_id[i] <- cl$clone_id[cl$cell_id == cndf$cell_id[i]]
}
cndf <- cndf[cndf$clone_id != "",]

## 1. Assign gene of interest
cndf$genei <- cndf[, j]

## 2. Set up clone ID factor levels
unique_clones <- sort(unique(as.character(cndf$clone_id)))
if (!exists("clone_levels")){
  clone_levels <- if ("0" %in% unique_clones) {
    c(setdiff(unique_clones, "0"), "0")
  } else {
    unique_clones
  }
}
cndf$clone_id <- factor(cndf$clone_id, levels = clone_levels)

## 3. Definition of axis range
axis_max <- ceiling(max(cndf$genei, na.rm = TRUE)/50)*50
breaks <- seq(0, axis_max, by = ceiling(axis_max/50))

## 4. Generate histogram counts for each clone
data_list <- split(cndf$genei, cndf$clone_id)
hist_list <- lapply(data_list, function(x) hist(x, breaks = breaks, plot = FALSE))
counts_matrix <- do.call(rbind, lapply(hist_list, function(h) h$counts))

## 5. Set bar colors (gray last if clone 0 present)
base_colors <- c("#f9d14a", "#00688B", "#ab3329", "#191970", "#88a0dc", "#e78429", "#7c4b73", "#ed968c")
bar_colors <- if ("0" %in% clone_levels) {
  c(base_colors[1:(length(data_list) - 1)], "gray")
} else {
  base_colors[1:length(data_list)]
}
rm(clone_levels)

## 6. Compute simulated bar centers to get xlim
bin_width <- breaks[2] - breaks[1]
n_bins <- ncol(counts_matrix)
bar_centers_sim <- seq(from = bin_width / 2, by = bin_width, length.out = n_bins)
right_edges <- bar_centers_sim + bin_width / 2
left_edge <- min(left_edges)
right_edge <- max(left_edges + bin_width)

## 7. Draw stacked barplot
pdf(paste0("Fig2h.left.", w, "_", j, ".pdf"), height=4, width=4)
bar_centers <- barplot(counts_matrix, beside = FALSE, col = bar_colors, border = NA, space = 0, xlim = c(left_edge, right_edge), xlab = "Per-cell copy number", ylab = "Number of cells", main = paste0(w, " (n = ", nrow(cndf), ") ", j, " CN by clone"), las = 1, xaxt = "n")

## 8. Draw X axis ticks at ~6 clean values
tick_step <- if (axis_max <= 50) 5 else 10
valid_ticks <- breaks[breaks %% tick_step == 0]
tick_labels <- round(valid_ticks, 1)
tick_idx <- match(tick_labels, breaks)
tick_positions <- c(0, right_edges)[tick_idx]

axis(
  side = 1,
  at = tick_positions,
  labels = tick_labels,
  cex.axis = 1
)

## 9. Legend
legend("topright", legend = paste0("clone ", names(data_list)), fill = bar_colors, bty = "n", border = NA)
dev.off()

# Panel H, middle
w <- "SPECTRUM-OV-052"
clonelist <- c("C", "A", "B")
chrs <- c("11", "13", "12")

qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)

df_final <- fread(paste0("./normalization/cn_dataframe/", w, ".normalized.gc_corrected.filters.annotated.csv"))
df_final <- df_final[
  normal_reads != 0 &
    segment_length >= 1000 &
    centroacro == "no" &
    blacklist == "no" &
    xfilter == "no"
]
df_final$cell_id <- str_sub(df_final$cell_id, -15)
df_final <- df_final[cell_id %in% qc$cell_id]

cl <- read.csv(paste0("./SPECTRUM-OV-052/SPECTRUM-OV-052_clusters_0.3_1000.csv"))
cl$cell_id <- str_sub(cl$cell_id, -15)
colnames(cl)[2] <- "clone_id"
cl$clone_id[cl$clone_id %in% c("C", "D", "E", "F", "H", "I")] <- "0"
cl$clone_id[cl$clone_id == "G"] <- "C"

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

isv <- read.csv(paste0("./hmf/gridss_", ifelse(indi$normal_wgs[indi$individual == w] == "yes", "somatic", "tumoronly"), "/", w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.annotated.bedpe"), header=T, as.is=T, sep="\t")
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

pdf(paste0("Fig2h.middle.", w, "_", paste(chrs, collapse = "-"), ".clones.", paste(clonelist, collapse = "-"), ".pdf"), width=6, height=4)
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

k <- 5
for (l in clonelist){
  df <- df_final[cell_id %in% cl$cell_id[cl$clone_id == l] & chr %in% chrs]
  cnvdf <- consensuscopynumber(df)
  if (k < max(cnvdf$copy, na.rm = T)){
    k <- max(cnvdf$copy, na.rm = T)
  }
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

for (l in clonelist){
  df <- df_final[cell_id %in% cl$cell_id[cl$clone_id == l] & chr %in% chrs]
  cnvdf <- consensuscopynumber(df)
  colnames(cnvdf) <- c("chromosome", "start", "end", "cell_id", "state", "totalCopyNumber")
  
  for (j in 1:nrow(cnvdf)){
    cnvdf$start[j] <- cnvdf$start[j] + deltadf$delta[deltadf$chrs == cnvdf$chromosome[j]]
    cnvdf$end[j] <- cnvdf$end[j] + deltadf$delta[deltadf$chrs == cnvdf$chromosome[j]]
  }
  cnvdf <- cnvdf[order(cnvdf$start, decreasing = F),]
  cnvdf$seglen <- cnvdf$end - cnvdf$start
  cnvdf <- cnvdf[cnvdf$seglen >= 500,]
  
  segments(cnvdf$start[cnvdf$start != cnvdf$end], cnvdf$totalCopyNumber[cnvdf$start != cnvdf$end], cnvdf$end[cnvdf$start != cnvdf$end], cnvdf$totalCopyNumber[cnvdf$start != cnvdf$end], col=base_colors[which(l == sort(clonelist))], lwd=1, lend=0)
  for (j in 2:nrow(cnvdf)){
    segments(cnvdf$start[j], cnvdf$totalCopyNumber[j-1], cnvdf$start[j], cnvdf$totalCopyNumber[j], col=base_colors[which(l == sort(clonelist))], lwd=1, lend=0)
  }
}

legend("topright", legend = paste0("clone ", sort(clonelist)), col = base_colors[1:length(sort(clonelist))], lty = 1, bty = "n", border = NA)

for (j in c("CCND1", "NCOA2")){
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

# Panel H, right: conceptual illustration

# Panel I, left -- we renamed subclones in the manuscript to follow the flow.
w <- "SPECTRUM-OV-083"
j <- "NCOA2"
clonelist <- c("G", "B", "H", "I", "D", "0")

## Loading the data
cndf <- read.csv(paste0("./normalization/oncogene_cn_final/", w, ".sv-seg-cn.oncogene.cn.csv"))
cndf$cell_id <- str_sub(cndf$cell_id, -15)
metrics <- read.csv(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
metrics$cell_id <- str_sub(metrics$cell_id, -15)
cndf <- cndf[cndf$cell_id %in% metrics$cell_id[metrics$decision == "include"],]

cl <- read.csv(paste0("./recluster_nov_24/tables/", w, "/", w, "_clusters_0.3_1000.csv"))
cl$cell_id <- str_sub(cl$cell_id, -15)
colnames(cl)[2] <- "clone_id"
cl$clone_id[cl$clone_id %in% c("A", "C", "E", "F")] <- "0"
cndf$clone_id <- ""
for (i in which(cndf$cell_id %in% cl$cell_id)){
  cndf$clone_id[i] <- cl$clone_id[cl$cell_id == cndf$cell_id[i]]
}
cndf <- cndf[cndf$clone_id != "",]

## 1. Assign gene of interest
cndf$genei <- cndf[, j]

## 2. Set up clone ID factor levels
unique_clones <- sort(unique(as.character(cndf$clone_id)))
if (!exists("clone_levels")){
  clone_levels <- if ("0" %in% unique_clones) {
    c(setdiff(unique_clones, "0"), "0")
  } else {
    unique_clones
  }
}
cndf$clone_id <- factor(cndf$clone_id, levels = clone_levels)

## 3. Definition of axis range
axis_max <- ceiling(max(cndf$genei, na.rm = TRUE)/50)*50
breaks <- seq(0, axis_max, by = ceiling(axis_max/50))

## 4. Generate histogram counts for each clone
data_list <- split(cndf$genei, cndf$clone_id)
hist_list <- lapply(data_list, function(x) hist(x, breaks = breaks, plot = FALSE))
counts_matrix <- do.call(rbind, lapply(hist_list, function(h) h$counts))

## 5. Set bar colors (gray last if clone 0 present)
base_colors <- c("#ab3329", "#f9d14a", "#00688B", "#191970", "#88a0dc", "#e78429", "#7c4b73", "#ed968c")
bar_colors <- if ("0" %in% clone_levels) {
  c(base_colors[1:(length(data_list) - 1)], "gray")
} else {
  base_colors[1:length(data_list)]
}
rm(clone_levels)

## 6. Compute simulated bar centers to get xlim
bin_width <- breaks[2] - breaks[1]
n_bins <- ncol(counts_matrix)
bar_centers_sim <- seq(from = bin_width / 2, by = bin_width, length.out = n_bins)
right_edges <- bar_centers_sim + bin_width / 2
left_edge <- min(left_edges)
right_edge <- max(left_edges + bin_width)

## 7. Draw stacked barplot
pdf(paste0("Fig2i.left.", w, "_", j, ".pdf"), height=4, width=4)
bar_centers <- barplot(counts_matrix, beside = FALSE, col = bar_colors, border = NA, space = 0, xlim = c(left_edge, right_edge), xlab = "Per-cell copy number", ylab = "Number of cells", main = paste0(w, " (n = ", nrow(cndf), ") ", j, " CN by clone"), las = 1, xaxt = "n")

## 8. Draw X axis ticks at ~6 clean values
tick_step <- if (axis_max <= 50) 5 else 10
valid_ticks <- breaks[breaks %% tick_step == 0]
tick_labels <- round(valid_ticks, 1)
tick_idx <- match(tick_labels, breaks)
tick_positions <- c(0, right_edges)[tick_idx]

axis(
  side = 1,
  at = tick_positions,
  labels = tick_labels,
  cex.axis = 1
)

## 9. Legend
legend("topright", legend = paste0("clone ", names(data_list)), fill = bar_colors, bty = "n", border = NA)
dev.off()

# Panel I, middle
w <- "SPECTRUM-OV-083"
clonelist <- c("G", "B", "H", "I", "D")
chrs <- c("8")

qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)

df_final <- fread(paste0("./normalization/cn_dataframe/", w, ".normalized.gc_corrected.filters.annotated.csv"))
df_final <- df_final[
  normal_reads != 0 &
    segment_length >= 1000 &
    centroacro == "no" &
    blacklist == "no" &
    xfilter == "no"
]
df_final$cell_id <- str_sub(df_final$cell_id, -15)
df_final <- df_final[cell_id %in% qc$cell_id]

cl <- read.csv("./recluster_nov_24/tables/SPECTRUM-OV-083/SPECTRUM-OV-083_clusters_0.3_1000.csv")
cl$cell_id <- str_sub(cl$cell_id, -15)
colnames(cl)[2] <- "clone_id"
cl$clone_id[cl$clone_id %in% c("E", "F")] <- "0"

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

isv <- read.csv(paste0("./hmf/gridss_", ifelse(indi$normal_wgs[indi$individual == w] == "yes", "somatic", "tumoronly"), "/", w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.annotated.bedpe"), header=T, as.is=T, sep="\t")
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

pdf(paste0("Fig2i.middle.", w, "_", paste(chrs, collapse = "-"), ".clones.", paste(clonelist, collapse = "-"), ".pdf"), width=6, height=4)
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

k <- 5
for (l in clonelist){
  df <- df_final[cell_id %in% cl$cell_id[cl$clone_id == l] & chr %in% chrs]
  cnvdf <- consensuscopynumber(df)
  if (k < max(cnvdf$copy, na.rm = T)){
    k <- max(cnvdf$copy, na.rm = T)
  }
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

for (l in clonelist){
  df <- df_final[cell_id %in% cl$cell_id[cl$clone_id == l] & chr %in% chrs]
  cnvdf <- consensuscopynumber(df)
  colnames(cnvdf) <- c("chromosome", "start", "end", "cell_id", "state", "totalCopyNumber")
  
  for (j in 1:nrow(cnvdf)){
    cnvdf$start[j] <- cnvdf$start[j] + deltadf$delta[deltadf$chrs == cnvdf$chromosome[j]]
    cnvdf$end[j] <- cnvdf$end[j] + deltadf$delta[deltadf$chrs == cnvdf$chromosome[j]]
  }
  cnvdf <- cnvdf[order(cnvdf$start, decreasing = F),]
  cnvdf$seglen <- cnvdf$end - cnvdf$start
  cnvdf <- cnvdf[cnvdf$seglen >= 500,]
  
  segments(cnvdf$start[cnvdf$start != cnvdf$end], cnvdf$totalCopyNumber[cnvdf$start != cnvdf$end], cnvdf$end[cnvdf$start != cnvdf$end], cnvdf$totalCopyNumber[cnvdf$start != cnvdf$end], col=base_colors[which(l == sort(clonelist))], lwd=1, lend=0)
  for (j in 2:nrow(cnvdf)){
    segments(cnvdf$start[j], cnvdf$totalCopyNumber[j-1], cnvdf$start[j], cnvdf$totalCopyNumber[j], col=base_colors[which(l == sort(clonelist))], lwd=1, lend=0)
  }
}

legend("topright", legend = paste0("clone ", sort(clonelist)), col = base_colors[1:length(sort(clonelist))], lty = 1, bty = "n", border = NA)

for (j in c("CCND1", "NCOA2")){
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

# Panel I, right: conceptual illustration

