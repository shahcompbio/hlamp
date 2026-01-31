# Figure 5
library(stringr)
library(psych)
indi <- read.csv("./metadata/dlp.cohort.summary.individual_111925.tsv", header=T, as.is=T, sep="\t")
east <- read.csv("./modeling/amplified_regions_summary_table_011026.csv")

# Panel A
aadf <- read.csv("./aa/aa.output.summary.123125.tsv", header=T, as.is=T, sep="\t")
df <- east[east$mean_cn > 6,]
df$aaac <- "not_available"
for (w in unique(df$individual)){
  tempdf <- as.data.frame(unique(str_sub(strsplit(paste(str_sub(aadf$Location[aadf$Individual == w], 2, -2), collapse = ", "), ", ", fixed = T)[[1]], 2, -2)))
  if (nrow(tempdf) != 0){
    colnames(tempdf) <- "reg"
    tempdf$chr <- ""
    tempdf$start <- NA
    tempdf$end <- NA
    tempdf$key <- NA
    for (i in 1:nrow(tempdf)){
      tempdf$chr[i] <- strsplit(tempdf$reg[i], ":", fixed = T)[[1]][1]
      tempdf$start[i] <- as.numeric(strsplit(strsplit(tempdf$reg[i], ":", fixed = T)[[1]][2], "-", fixed = T)[[1]][1])
      tempdf$end[i] <- as.numeric(strsplit(strsplit(tempdf$reg[i], ":", fixed = T)[[1]][2], "-", fixed = T)[[1]][2])
      tempdf$key[i] <- paste(which(grepl(tempdf$reg[i], aadf$Location[aadf$Individual == w & !is.na(aadf$Classification)])), collapse = ";")
    }
    for (i in which(df$individual == w)){
      if (nrow(tempdf[tempdf$chr == df$chromosome[i] & tempdf$start <= df$end[i] & tempdf$end >= df$start[i],]) != 0){
        if ("ecDNA" %in% aadf$Classification[aadf$Individual == w][as.numeric(strsplit(tempdf$key[tempdf$chr == df$chromosome[i] & tempdf$start <= df$end[i] & tempdf$end >= df$start[i]], ";", fixed = T)[[1]])]){
          df$aaac[i] <- "ecDNA"
        } else {
          df$aaac[i] <- "ICamp"
        }
      }
    }
  }
}
df$aaac <- factor(df$aaac, levels = c("ecDNA", "ICamp", "not_available"))

## 2x2 table values
length(unique(df$individual)) # n=63 with HLAMPs
length(unique(df$individual[df$decision == "ecDNA"])) # 11
length(unique(df$individual[df$decision == "ecDNA" & df$aaac == "ecDNA"])) # 8 (+1 -- 2765_2 in mouse run) --> 9
length(unique(aadf$Individual[aadf$Classification == "ecDNA"])) # 40 -- 40 cases with ecDNA positive by AA
nrow(indi) - length(unique(aadf$Individual[aadf$Classification == "ecDNA"])) - 2

## Cohen's kappa
cohen.kappa(matrix(c(9,31,2,60), nrow = 2))

# Panel B, left
mat <- matrix(NA, nrow = 7, ncol = 3)
row.names(mat) <- rev(names(table(df$cancer_type)))
colnames(mat) <- c("ecDNA", "ICamp", "not_available")

for (i in 1:nrow(mat)){
  if (rownames(mat)[i] == "HGSOC"){
    mat["HGSOC",] <- unname(table(df$aaac[df$cancer_type == "HGSOC" & df$decision == "ICamp"]))
  } else {
    mat[i,] <- unname(table(df$aaac[df$cancer_type == rownames(mat)[i] & df$ecdna == "no"]))
  }
}
pdf("Fig5b.left.pdf", height = 4, width = 3)
barplot(t(mat), horiz = T, las=1, xlab = "# of ICamp regions", xlim = c(0,260))
dev.off()

tempdf <- as.data.frame(mat)
mat <- matrix(NA, nrow = 7, ncol = 3)
row.names(mat) <- rev(names(table(df$cancer_type)))
colnames(mat) <- c("ecDNA", "ICamp", "not_available")

for (i in 1:nrow(mat)){
  if (rownames(mat)[i] == "HGSOC"){
    mat["HGSOC",] <- unname(table(df$aaac[df$cancer_type == "HGSOC" & df$decision == "ecDNA"]))
  } else {
    mat[i,] <- unname(table(df$aaac[df$cancer_type == rownames(mat)[i] & df$ecdna == "yes"]))
  }
}
pdf("Fig5b.right.pdf", height = 4, width = 3)
barplot(t(mat), horiz = T, las=1, xlab = "# of ecDNA regions", xlim = c(0,260))
dev.off()

# Panel C
df <- east[east$mean_cn >= 6 & !is.na(east$rightcn_s) & !is.na(east$leftcn_s) & !is.na(east$rightcn_e) & !is.na(east$leftcn_e),]
df$decision <- factor(df$decision, levels = c("ecDNA", "ICamp"))
df$border_cn_gap <- ((df$rightcn_s - df$leftcn_s) + (df$leftcn_e - df$rightcn_e))/2

pdf("Fig5c.pdf", height=4, width=2.7)
options(scipen = 0)
boxplot(df$border_cn_gap ~ df$decision, frame = F, las = 1, ylim = c(0,200), xlab = "", ylab = "Mean copy number gap at the HLAMP border", outline = F, col = c("#F0E68C", "#FFBBFF"), main = paste0("p = ", t.test(df$border_cn_gap ~ df$decision)$p.value), lwd = 2/3)
stripchart(df$border_cn_gap[!(df$clonal90 == "no" & df$decision == "ecDNA") & !(df$individual %in% c("NCI-H82", "COLO320-HSR", "GBM39-HSR"))] ~ df$decision[!(df$clonal90 == "no" & df$decision == "ecDNA") & !(df$individual %in% c("NCI-H82", "COLO320-HSR", "GBM39-HSR"))], pch = 19, col = rgb(0,0,0,.1), method = "jitter", vertical = T, add = T)
stripchart(df$border_cn_gap[df$clonal90 == "no" & df$decision == "ecDNA"] ~ df$decision[df$clonal90 == "no" & df$decision == "ecDNA"], pch = 20, col = "#FF9912", method = "jitter", vertical = T, add = T)
stripchart(df$border_cn_gap[df$individual %in% c("NCI-H82", "COLO320-HSR", "GBM39-HSR")] ~ df$decision[df$individual %in% c("NCI-H82", "COLO320-HSR", "GBM39-HSR")], pch = 20, col = "#8B008B", method = "jitter", vertical = T, add = T)
legend("topright", legend = c("Subclonal ecDNA", "HSR"), pch = 20, col = c("#FF9912", "#8B008B"), border = F)
options(scipen = 999)
dev.off()

# Panel D, left
w <- "COLO320DM"
j <- "MYC"
cndf <- read.csv(paste0("./normalization/oncogene_cn_final/", w, ".sv-seg-cn.oncogene.cn.csv"))
cndf$cell_id <- str_sub(cndf$cell_id, -15)
qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)
cndf <- cndf[cndf$cell_id %in% qc$cell_id[qc$decision == "include"],]

infoline <- cndf[,j]
infoline <- infoline[!is.na(infoline) & infoline >0]

pdf(paste0("Fig5d.left.", w, ".", j, ".pdf"), height=4, width=3.5)
hist(infoline, las=1, breaks = seq(0,1100,1100/50), ylim = c(0,0.036), xlab = "Per cell copy number", main = paste0(w, " (n=", length(infoline), ") ", j), col = "#FF82AB", border = F, prob = T)
lines(density(infoline), lwd=2, col="#B0171F")
dev.off()

w <- "COLO320-HSR"
j <- "MYC"
cndf <- read.csv(paste0("./normalization/oncogene_cn_final/", w, ".sv-seg-cn.oncogene.cn.csv"))
cndf$cell_id <- str_sub(cndf$cell_id, -15)
qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)
cndf <- cndf[cndf$cell_id %in% qc$cell_id[qc$decision == "include"],]

infoline <- cndf[,j]
infoline <- infoline[!is.na(infoline) & infoline >0]

pdf(paste0("Fig5d.left.", w, ".", j, ".pdf"), height=4, width=3.5)
hist(infoline, las=1, breaks = seq(0,1100,1100/50), ylim = c(0,0.036), xlab = "Per cell copy number", main = paste0(w, " (n=", length(infoline), ") ", j), col = "#1E90FF", border = F, prob = T)
lines(density(infoline), lwd=2, col="#0000FF")
dev.off()

# Panel D, right
w <- "GBM39-DM"
j <- "EGFR"
cndf <- read.csv(paste0("./normalization/oncogene_cn_final/", w, ".sv-seg-cn.oncogene.cn.csv"))
cndf$cell_id <- str_sub(cndf$cell_id, -15)
qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)
cndf <- cndf[cndf$cell_id %in% qc$cell_id[qc$decision == "include"],]

infoline <- cndf[,j]
infoline <- infoline[!is.na(infoline) & infoline >0]

pdf(paste0("Fig5d.right.", w, ".", j, ".pdf"), height=4, width=3.5)
hist(infoline, las=1, breaks = seq(0,700,700/50), ylim = c(0,0.032), xlab = "Per cell copy number", main = paste0(w, " (n=", length(infoline), ") ", j), col = "#FF82AB", border = F, prob = T)
lines(density(infoline), lwd=2, col="#B0171F")
dev.off()

w <- "GBM39-HSR"
j <- "EGFR"
cndf <- read.csv(paste0("./normalization/oncogene_cn_final/", w, ".sv-seg-cn.oncogene.cn.csv"))
cndf$cell_id <- str_sub(cndf$cell_id, -15)
qc <- fread(paste0("./metadata/doublet_invariant_filtered/", w, "_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)
cndf <- cndf[cndf$cell_id %in% qc$cell_id[qc$decision == "include"],]

infoline <- cndf[,j]
infoline <- infoline[!is.na(infoline) & infoline >0]

pdf(paste0("Fig5d.right.", w, ".", j, ".pdf"), height=4, width=3.5)
hist(infoline, las=1, breaks = seq(0,700,700/50), ylim = c(0,0.032), xlab = "Per cell copy number", main = paste0(w, " (n=", length(infoline), ") ", j), col = "#1E90FF", border = F, prob = T)
lines(density(infoline), lwd=2, col="#0000FF")
dev.off()

# Panel E: FISH images

# Panel F
ptlist <- c("COLO320DM", "COLO320-HSR")
k <- 300
regions <- "6-1-1000000;8-127000001-131000000;13-28000001-29000000"

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

pdf("Fig5f.pdf", height=4, width=7)
plot(0, type='n', xlim=c(0, df$add_value[nrow(df)] + df$seglen[nrow(df)]), ylim=c(0,1.1*k), xaxt = 'n', frame=F, las=1, xlab=paste0("Position on chromosome ", paste(unique(df$chr), collapse = "-"), " in Mbp"), ylab="Total copy number", main = "COLO320")
for (j in 1:nrow(df)){
  axis(1, at = c(df$add_value[j], df$add_value[j] + df$seglen[j]), tick = c(df$add_value[j]+1, df$add_value[j] + df$seglen[j]), labels = c(df$start[j], df$end[j]), las=2)
}

isv <- read.csv("./sc_sv_genotyping/COLO320_merged.annotated.bedpe", header=T, as.is=T, sep="\t")

svdf <- isv
setDT(svdf)
svdf[, chr1 := as.character(chr1)]
svdf[, chr2 := as.character(chr2)]
swap <- chr_rank(svdf$chr2) < chr_rank(svdf$chr1)
svdf[swap, c("chr1","start1","end1","ori1","chr2","start2","end2","ori2") :=
       .(chr2, start2, end2, ori2, chr1, start1, end1, ori1)]

svdf$pos1 <- round((svdf$start1 + svdf$end1)/2)
svdf$pos2 <- round((svdf$start2 + svdf$end2)/2)
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

svdf[, c("reg1", "reg2") := lapply(.SD, as.character), .SDcols = c("reg1", "reg2")]
df[, regs := as.character(regs)]

svdf <- as.data.frame(svdf)
svdf$plotpos1 <- svdf$pos1
svdf$plotpos2 <- svdf$pos2
for (j in 1:nrow(df)){
  svdf$plotpos1[svdf$reg1 == df$regs[j] & !is.na(svdf$reg1)] <- svdf$plotpos1[svdf$reg1 == df$regs[j] & !is.na(svdf$reg1)] - df$start[j] + df$add_value[j]
  svdf$plotpos2[svdf$reg2 == df$regs[j] & !is.na(svdf$reg2)] <- svdf$plotpos2[svdf$reg2 == df$regs[j] & !is.na(svdf$reg2)] - df$start[j] + df$add_value[j]
}
setDT(svdf)

## SV plotting
theta=seq(0,pi, len=100)
chrs <- df$chr
svdf <- as.data.frame(svdf)
svdf <- svdf[svdf$svname != "hsrsv_4",] # This SV is a pericentromeric artifact with extreme number of read support -- was manually excluded
svdf <- svdf[!(svdf$dm_alt_reads == 0 & svdf$hsr_alt_reads == 0),]

svdf$color <- ""
svdf$color[svdf$dm_alt_reads >= 10 & svdf$hsr_alt_reads == 0] <- "#B0171F" # DM-specific SVs -- SVs that occurred after the clonal divergence, specifically in the DM cell line
svdf$color[svdf$dm_alt_reads == 0 & svdf$hsr_alt_reads >= 10] <- "#0000FF" # HSR-specific SVs -- SVs that occurred after the clonal divergence, specifically in the HSR cell line
svdf$color[svdf$color == ""] <- "#FFB90F" # common SVs -- SVs that happened before the clonal divergence
svdf$composite_vf <- 0
svdf$composite_vf[svdf$dm_alt_reads >= 10 & svdf$hsr_alt_reads == 0] <- svdf$dm_alt_reads[svdf$dm_alt_reads >= 10 & svdf$hsr_alt_reads == 0]
svdf$composite_vf[svdf$dm_alt_reads == 0 & svdf$hsr_alt_reads >= 10] <- svdf$hsr_alt_reads[svdf$dm_alt_reads == 0 & svdf$hsr_alt_reads >= 10]
svdf$composite_vf[svdf$color == "#FFB90F"] <- (svdf$dm_alt_reads[svdf$color == "#FFB90F"] + svdf$hsr_alt_reads[svdf$color == "#FFB90F"])/2
svdf$composite_vf <- log(svdf$composite_vf)

isvin <- svdf[svdf$chr1 %in% chrs & svdf$chr2 %in% chrs,]
isvin$rad=abs(isvin$plotpos2-isvin$plotpos1)/2
isvout <- svdf[!(svdf$chr1 %in% chrs) | !(svdf$chr2 %in% chrs),]

max_val <- max(svdf$composite_vf)

par(new = T)
plot(0, type='n', xlim=c(0, df$add_value[nrow(df)] + df$seglen[nrow(df)]), ylim=c(0,1.1*max_val), frame=F, axes = F, xlab = "", ylab = "")
axis(4, las=1)
mtext("Log(read support)", side=4, line=3)

for (j in 1:nrow(isvin)){
  if (isvin$strand1[j] == "+" & isvin$strand2[j] == "-"){
    x <- isvin$rad[j]*cos(theta)+((isvin$plotpos1[j]+isvin$plotpos2[j])/2)
    y <- max_val*0.1*sin(theta)+isvin$composite_vf[j]
  } else if (isvin$strand1[j] == "-" & isvin$strand2[j] == "+"){
    x <- isvin$rad[j]*cos(theta)+((isvin$plotpos1[j]+isvin$plotpos2[j])/2)
    y <- -max_val*0.1*sin(theta)+isvin$composite_vf[j]
  } else if (isvin$strand1[j] == "-" & isvin$strand2[j] == "-"){
    x <- isvin$rad[j]*cos(theta)+((isvin$plotpos1[j]+isvin$plotpos2[j])/2)
    y <- max_val*0.1*sin(theta)+isvin$composite_vf[j]
  } else if (isvin$strand1[j] == "+" & isvin$strand2[j] == "+"){
    x <- isvin$rad[j]*cos(theta)+((isvin$plotpos1[j]+isvin$plotpos2[j])/2)
    y <- -max_val*0.1*sin(theta)+isvin$composite_vf[j]
  }
  lines(x,y,col=isvin$color[j], xpd=TRUE, lwd=1/3)
  segments(isvin$plotpos1[j], 0, isvin$plotpos1[j], isvin$composite_vf[j], lty=1, xpd=TRUE, col=isvin$color[j], lwd=1/3)
  segments(isvin$plotpos2[j], 0, isvin$plotpos2[j], isvin$composite_vf[j], lty=1, xpd=TRUE, col=isvin$color[j], lwd=1/3)
}

isvout$composite_vf[isvout$composite_vf < 1] <- 1
for (j in 1:nrow(isvout)){
  if (isvout$chr1[j] %in% chrs){
    arrows(isvout$plotpos1[j], 0, isvout$plotpos1[j], isvout$composite_vf[j], lty=1, xpd=TRUE, col=isvout$color[j], code=2, length = 0.1, lwd=1/3)
  } else {
    arrows(isvout$plotpos2[j], 0, isvout$plotpos2[j], isvout$composite_vf[j], lty=1, xpd=TRUE, col=isvout$color[j], code=2, length = 0.1, lwd=1/3)
  }
}

## CNA -- COLO320DM
qc <- fread(paste0("./metadata/doublet_invariant_filtered/COLO320DM_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)

df_final <- fread(paste0("./normalization/cn_dataframe/COLO320DM.normalized.gc_corrected.filters.annotated.csv"))
df_final$cell_id <- str_sub(df_final$cell_id, -15)
df_final <- df_final[
  normal_reads != 0 &
    segment_length >= 1000 &
    centroacro == "no" &
    blacklist == "no" &
    xfilter == "no" &
    cell_id %in% qc$cell_id[qc$decision == "include"]
]

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
    .rowid = x..rowid # Keep the original row ID
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

cnvdf <- consensuscopynumber(as.data.frame(df_final[,c("chr", "start", "end", "reads", "copy", "state", "cell_id")]))

cnvdf$plot.start <- NA
cnvdf$plot.end <- NA
for (j in unique(cnvdf$chr)){
  cnvdf$plot.start[cnvdf$chr == j] <- cnvdf$start[cnvdf$chr == j] - df$start[df$chr == j] + df$add_value[df$chr == j]
  cnvdf$plot.end[cnvdf$chr == j] <- cnvdf$end[cnvdf$chr == j] - df$start[df$chr == j] + df$add_value[df$chr == j]
}

cnvdf <- cnvdf[,c("chr", "plot.start", "plot.end", "copy", "state")]
colnames(cnvdf) <- c("chromosome", "plotstart", "plotend", "CopyNumber", "state")
cnvdf$state[cnvdf$state >= 11] <- 11
cnvdf$state[cnvdf$state < 0] <- 0
cnvdf$plotstart <- as.numeric(cnvdf$plotstart)
cnvdf$plotend <- as.numeric(cnvdf$plotend)

par(new = T)
plot(0, type='n', xlim=c(0, df$add_value[nrow(df)] + df$seglen[nrow(df)]), ylim=c(0,1.1*k), frame=F, axes = F, xlab = "", ylab = "")
rect(xleft = cnvdf$plotstart, ybottom = 0, xright = cnvdf$plotend, ytop = cnvdf$CopyNumber, col="#B0171F55", border = NA)

## CNA -- COLO320-HSR
qc <- fread(paste0("./metadata/doublet_invariant_filtered/COLO320-HSR_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)

df_final <- fread(paste0("./normalization/cn_dataframe/COLO320-HSR.normalized.gc_corrected.filters.annotated.csv"))
df_final$cell_id <- str_sub(df_final$cell_id, -15)
df_final <- df_final[
  normal_reads != 0 &
    segment_length >= 1000 &
    centroacro == "no" &
    blacklist == "no" &
    xfilter == "no" &
    cell_id %in% qc$cell_id[qc$decision == "include"]
]

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
  v[is.na(v)] <- 9999 
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

cnvdf <- consensuscopynumber(as.data.frame(df_final[,c("chr", "start", "end", "reads", "copy", "state", "cell_id")]))

cnvdf$plot.start <- NA
cnvdf$plot.end <- NA
for (j in unique(cnvdf$chr)){
  cnvdf$plot.start[cnvdf$chr == j] <- cnvdf$start[cnvdf$chr == j] - df$start[df$chr == j] + df$add_value[df$chr == j]
  cnvdf$plot.end[cnvdf$chr == j] <- cnvdf$end[cnvdf$chr == j] - df$start[df$chr == j] + df$add_value[df$chr == j]
}

cnvdf <- cnvdf[,c("chr", "plot.start", "plot.end", "copy", "state")]
colnames(cnvdf) <- c("chromosome", "plotstart", "plotend", "CopyNumber", "state")
cnvdf$state[cnvdf$state >= 11] <- 11
cnvdf$state[cnvdf$state < 0] <- 0
cnvdf$plotstart <- as.numeric(cnvdf$plotstart)
cnvdf$plotend <- as.numeric(cnvdf$plotend)

rect(xleft = cnvdf$plotstart, ybottom = 0, xright = cnvdf$plotend, ytop = cnvdf$CopyNumber, col="#0000FF55", border = NA)
dev.off()

# Panel G
ptlist <- c("GBM39-DM", "GBM39-HSR")
k <- 200

regions <- "7-50000001-105000000;9-90000001-141213431"
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

pdf("Fig5g.pdf", height=4, width=7)
plot(0, type='n', xlim=c(0, df$add_value[nrow(df)] + df$seglen[nrow(df)]), ylim=c(0,1.1*k), xaxt = 'n', frame=F, las=1, xlab=paste0("Position on chromosome ", paste(unique(df$chr), collapse = "-"), " in Mbp"), ylab="Total copy number", main = "GBM39")
for (j in 1:nrow(df)){
  axis(1, at = c(df$add_value[j], df$add_value[j] + df$seglen[j]), tick = c(df$add_value[j]+1, df$add_value[j] + df$seglen[j]), labels = c(df$start[j], df$end[j]), las=2)
}

isv <- read.csv("./sc_sv_genotyping/GBM39_merged.annotated.bedpe", header=T, as.is=T, sep="\t")

svdf <- isv
setDT(svdf)
svdf[, chr1 := as.character(chr1)]
svdf[, chr2 := as.character(chr2)]
swap <- chr_rank(svdf$chr2) < chr_rank(svdf$chr1)
svdf[swap, c("chr1","start1","end1","ori1","chr2","start2","end2","ori2") :=
       .(chr2, start2, end2, ori2, chr1, start1, end1, ori1)]

svdf$pos1 <- round((svdf$start1 + svdf$end1)/2)
svdf$pos2 <- round((svdf$start2 + svdf$end2)/2)
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

svdf[, c("reg1", "reg2") := lapply(.SD, as.character), .SDcols = c("reg1", "reg2")]
df[, regs := as.character(regs)]

svdf <- as.data.frame(svdf)
svdf$plotpos1 <- svdf$pos1
svdf$plotpos2 <- svdf$pos2
for (j in 1:nrow(df)){
  svdf$plotpos1[svdf$reg1 == df$regs[j] & !is.na(svdf$reg1)] <- svdf$plotpos1[svdf$reg1 == df$regs[j] & !is.na(svdf$reg1)] - df$start[j] + df$add_value[j]
  svdf$plotpos2[svdf$reg2 == df$regs[j] & !is.na(svdf$reg2)] <- svdf$plotpos2[svdf$reg2 == df$regs[j] & !is.na(svdf$reg2)] - df$start[j] + df$add_value[j]
}
setDT(svdf)

## SV plotting
theta=seq(0,pi, len=100)
chrs <- df$chr
svdf <- as.data.frame(svdf)
svdf <- svdf[!(svdf$dm_alt_reads == 0 & svdf$hsr_alt_reads == 0),]

svdf$color <- ""
svdf$color[svdf$dm_alt_reads >= 10 & svdf$hsr_alt_reads == 0] <- "#B0171F" # DM-specific SVs -- SVs that occurred after the clonal divergence, specifically in the DM cell line
svdf$color[svdf$dm_alt_reads == 0 & svdf$hsr_alt_reads >= 10] <- "#0000FF" # HSR-specific SVs -- SVs that occurred after the clonal divergence, specifically in the HSR cell line
svdf$color[svdf$color == ""] <- "#FFB90F" # common SVs -- SVs that happened before the clonal divergence
svdf$composite_vf <- 0
svdf$composite_vf[svdf$dm_alt_reads >= 10 & svdf$hsr_alt_reads == 0] <- svdf$dm_alt_reads[svdf$dm_alt_reads >= 10 & svdf$hsr_alt_reads == 0]
svdf$composite_vf[svdf$dm_alt_reads == 0 & svdf$hsr_alt_reads >= 10] <- svdf$hsr_alt_reads[svdf$dm_alt_reads == 0 & svdf$hsr_alt_reads >= 10]
svdf$composite_vf[svdf$color == "#FFB90F"] <- (svdf$dm_alt_reads[svdf$color == "#FFB90F"] + svdf$hsr_alt_reads[svdf$color == "#FFB90F"])/2
svdf$composite_vf <- log(svdf$composite_vf)

isvin <- svdf[svdf$chr1 %in% chrs & svdf$chr2 %in% chrs,]
isvin$rad=abs(isvin$plotpos2-isvin$plotpos1)/2
isvout <- svdf[!(svdf$chr1 %in% chrs) | !(svdf$chr2 %in% chrs),]

max_val <- max(svdf$composite_vf)

par(new = T)
plot(0, type='n', xlim=c(0, df$add_value[nrow(df)] + df$seglen[nrow(df)]), ylim=c(0,1.1*max_val), frame=F, axes = F, xlab = "", ylab = "")
axis(4, las=1)
mtext("Log(read support)", side=4, line=3)

for (j in 1:nrow(isvin)){
  if (isvin$strand1[j] == "+" & isvin$strand2[j] == "-"){
    x <- isvin$rad[j]*cos(theta)+((isvin$plotpos1[j]+isvin$plotpos2[j])/2)
    y <- max_val*0.1*sin(theta)+isvin$composite_vf[j]
  } else if (isvin$strand1[j] == "-" & isvin$strand2[j] == "+"){
    x <- isvin$rad[j]*cos(theta)+((isvin$plotpos1[j]+isvin$plotpos2[j])/2)
    y <- -max_val*0.1*sin(theta)+isvin$composite_vf[j]
  } else if (isvin$strand1[j] == "-" & isvin$strand2[j] == "-"){
    x <- isvin$rad[j]*cos(theta)+((isvin$plotpos1[j]+isvin$plotpos2[j])/2)
    y <- max_val*0.1*sin(theta)+isvin$composite_vf[j]
  } else if (isvin$strand1[j] == "+" & isvin$strand2[j] == "+"){
    x <- isvin$rad[j]*cos(theta)+((isvin$plotpos1[j]+isvin$plotpos2[j])/2)
    y <- -max_val*0.1*sin(theta)+isvin$composite_vf[j]
  }
  lines(x,y,col=isvin$color[j], xpd=TRUE, lwd=1/3)
  segments(isvin$plotpos1[j], 0, isvin$plotpos1[j], isvin$composite_vf[j], lty=1, xpd=TRUE, col=isvin$color[j], lwd=1/3)
  segments(isvin$plotpos2[j], 0, isvin$plotpos2[j], isvin$composite_vf[j], lty=1, xpd=TRUE, col=isvin$color[j], lwd=1/3)
}

isvout$composite_vf[isvout$composite_vf < 1] <- 1
for (j in 1:nrow(isvout)){
  if (isvout$chr1[j] %in% chrs){
    arrows(isvout$plotpos1[j], 0, isvout$plotpos1[j], isvout$composite_vf[j], lty=1, xpd=TRUE, col=isvout$color[j], code=2, length = 0.1, lwd=1/3)
  } else {
    arrows(isvout$plotpos2[j], 0, isvout$plotpos2[j], isvout$composite_vf[j], lty=1, xpd=TRUE, col=isvout$color[j], code=2, length = 0.1, lwd=1/3)
  }
}

## CNA -- GBM39-DM
qc <- fread(paste0("./metadata/doublet_invariant_filtered/GBM39-DM_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)

df_final <- fread(paste0("./normalization/cn_dataframe/GBM39-DM.normalized.gc_corrected.filters.annotated.csv"))
df_final$cell_id <- str_sub(df_final$cell_id, -15)
df_final <- df_final[
  normal_reads != 0 &
    segment_length >= 1000 &
    #    centroacro == "no" &
    blacklist == "no" &
    xfilter == "no" &
    cell_id %in% qc$cell_id[qc$decision == "include"]
]

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

cnvdf <- consensuscopynumber(as.data.frame(df_final[,c("chr", "start", "end", "reads", "copy", "state", "cell_id")]))

cnvdf$plot.start <- NA
cnvdf$plot.end <- NA
for (j in unique(cnvdf$chr)){
  cnvdf$plot.start[cnvdf$chr == j] <- cnvdf$start[cnvdf$chr == j] - df$start[df$chr == j] + df$add_value[df$chr == j]
  cnvdf$plot.end[cnvdf$chr == j] <- cnvdf$end[cnvdf$chr == j] - df$start[df$chr == j] + df$add_value[df$chr == j]
}

cnvdf <- cnvdf[,c("chr", "plot.start", "plot.end", "copy", "state")]
colnames(cnvdf) <- c("chromosome", "plotstart", "plotend", "CopyNumber", "state")
cnvdf$state[cnvdf$state >= 11] <- 11
cnvdf$state[cnvdf$state < 0] <- 0
cnvdf$plotstart <- as.numeric(cnvdf$plotstart)
cnvdf$plotend <- as.numeric(cnvdf$plotend)

par(new = T)
plot(0, type='n', xlim=c(0, df$add_value[nrow(df)] + df$seglen[nrow(df)]), ylim=c(0,1.1*k), xaxt = 'n', frame=F, xlab = "", ylab = "", axes = F)
rect(xleft = cnvdf$plotstart, ybottom = 0, xright = cnvdf$plotend, ytop = cnvdf$CopyNumber, col="#B0171F55", border = NA)

## CNA -- GBM39-HSR
qc <- fread(paste0("./metadata/doublet_invariant_filtered/GBM39-HSR_qc.filtered.csv"))
qc$cell_id <- str_sub(qc$cell_id, -15)

df_final <- fread(paste0("./normalization/cn_dataframe/GBM39-HSR.normalized.gc_corrected.filters.annotated.csv"))
df_final$cell_id <- str_sub(df_final$cell_id, -15)
df_final <- df_final[
  normal_reads != 0 &
    segment_length >= 1000 &
    #centroacro == "no" &
    blacklist == "no" &
    xfilter == "no" &
    cell_id %in% qc$cell_id[qc$decision == "include"]
]

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

cnvdf <- consensuscopynumber(as.data.frame(df_final[,c("chr", "start", "end", "reads", "copy", "state", "cell_id")]))

cnvdf$plot.start <- NA
cnvdf$plot.end <- NA
for (j in unique(cnvdf$chr)){
  cnvdf$plot.start[cnvdf$chr == j] <- cnvdf$start[cnvdf$chr == j] - df$start[df$chr == j] + df$add_value[df$chr == j]
  cnvdf$plot.end[cnvdf$chr == j] <- cnvdf$end[cnvdf$chr == j] - df$start[df$chr == j] + df$add_value[df$chr == j]
}

cnvdf <- cnvdf[,c("chr", "plot.start", "plot.end", "copy", "state")]
colnames(cnvdf) <- c("chromosome", "plotstart", "plotend", "CopyNumber", "state")
cnvdf$state[cnvdf$state >= 11] <- 11
cnvdf$state[cnvdf$state < 0] <- 0
cnvdf$plotstart <- as.numeric(cnvdf$plotstart)
cnvdf$plotend <- as.numeric(cnvdf$plotend)

rect(xleft = cnvdf$plotstart, ybottom = 0, xright = cnvdf$plotend, ytop = cnvdf$CopyNumber, col="#0000FF55", border = NA)

dev.off()




