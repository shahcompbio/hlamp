library(dplyr)
library(stringr)
library(signals)
library(data.table)
library(anndata)

args <- commandArgs(trailingOnly = TRUE)
w <- args

indi <- read.csv("./metadata/dlp.cohort.summary.individual.tsv", header=T, as.is=T, sep="\t")
genedf <- read.csv("./references/genedf.GRCh37.all.genes.txt", header = T, as.is=T, sep="\t")
drvamp <- read.csv("./resources/v5_34/ref/37/common/DriverGenePanel.37.tsv", header = T, as.is=T, sep="\t")
drvamp <- drvamp$gene[drvamp$reportAmplification == "true"]
drvamp <- drvamp[!(drvamp %in% c("TG", "RSF1"))]
drvamp <- c(drvamp, "PAK1")
drvamp <- genedf[genedf$gene %in% drvamp,]

hg19 <- read.csv("./references/hg19.coordinates.csv", header=T, as.is=T)
hg19$centro_mid <- round((hg19$centro_start+hg19$centro_end)/2)

cnvdf <- read.csv(paste0("./analysis/normalization/subclonal_cn_v4/", w, "_clusters_0.7_1000_subclonal_cn.gc_corrected.csv"))
cnvdf <- cnvdf[,colnames(cnvdf)[colnames(cnvdf) != "length"]]
segs <- read.csv(ifelse(w %in% indi$individual[indi$normal_wgs == "yes"], paste0("./analysis/segs/somatic/", w, "/segment.sv.tsv"), paste0("./analysis/segs/tumoronly/", w, "/segment.sv.tsv")), header=T, as.is=T, sep="\t")
cnvdf <- cnvdf %>% left_join(segs %>% select(chr, start, end, ssv, esv, svf, evf, seglen), by = c("chr", "start", "end"))

isv <- read.csv(ifelse(w %in% indi$individual[indi$normal_wgs == "yes"], paste0("./analysis/hmf/gridss_somatic/", w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.bedpe"), paste0("./analysis/hmf/gridss_tumoronly/", w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.bedpe")), header=F, as.is=T, sep="\t")
isv$type <- ""
isv$type[isv$V10 == "+" & isv$V11 == "-"] <- "DEL"
isv$type[isv$V10 == "-" & isv$V11 == "+"] <- "DUP"
isv$type[isv$V10 == "+" & isv$V11 == "+"] <- "t2tINV"
isv$type[isv$V10 == "-" & isv$V11 == "-"] <- "h2hINV"
isv$type[isv$V2 != isv$V5] <- "TRA"

cnvdf$ssv_type <- ""
for (i in which(cnvdf$ssv != "")){
  if (nrow(isv[isv$V8 == paste0(str_sub(cnvdf$ssv[i], 1, -2), "h") | isv$V8 == paste0(str_sub(cnvdf$ssv[i], 1, -2), "o"),]) != 0){
    cnvdf$ssv_type[i] <- unique(isv$type[isv$V8 == paste0(str_sub(cnvdf$ssv[i], 1, -2), "h") | isv$V8 == paste0(str_sub(cnvdf$ssv[i], 1, -2), "o")])
  }
}
cnvdf$esv_type <- ""
for (i in which(cnvdf$esv != "")){
  if (nrow(isv[isv$V8 == paste0(str_sub(cnvdf$esv[i], 1, -2), "h") | isv$V8 == paste0(str_sub(cnvdf$esv[i], 1, -2), "o"),]) != 0){
    cnvdf$esv_type[i] <- unique(isv$type[isv$V8 == paste0(str_sub(cnvdf$esv[i], 1, -2), "h") | isv$V8 == paste0(str_sub(cnvdf$esv[i], 1, -2), "o")])
  }
}
cnvdf$ssv_type[cnvdf$ssv_type == "" & cnvdf$ssv != "" & !(cnvdf$ssv %in% c("centromere", "telomere"))] <- "SBE"
cnvdf$esv_type[cnvdf$esv_type == "" & cnvdf$esv != "" & !(cnvdf$esv %in% c("centromere", "telomere"))] <- "SBE"

cnvdf$arm <- "no"
arm_idx <- "no"
for (i in 1:nrow(cnvdf)){
  if (cnvdf$end[i] > hg19$centro_mid[hg19$chr == cnvdf$chr[i]]){
    arm_idx <- "q"
  } else {
    arm_idx <- "p"
  }
  cnvdf$arm[i] <- paste0(cnvdf$chr[i], arm_idx)
}

cnvdf$cnbase <- 1
for (i in unique(cnvdf$arm)[!(unique(cnvdf$arm) %in% c("13p", "14p", "15p", "21p", "22p"))]){
  infodf1 <- NULL
  for (k in sort(unique(round(cnvdf$copy[cnvdf$arm == i])))){
    lenvalue <- sum(cnvdf$seglen[cnvdf$arm == i & round(cnvdf$copy) == k], na.rm = T)
    infoline1 <- cbind(i, k, lenvalue)
    infodf1 <- rbind(infodf1, infoline1)
  }
  infodf1 <- as.data.frame(infodf1)
  colnames(infodf1) <- c("arm", "total_cn", "lenvalue")
  infodf1$arm <- as.character(infodf1$arm)
  infodf1$total_cn <- as.numeric(as.character(infodf1$total_cn))
  infodf1$lenvalue <- as.numeric(as.character(infodf1$lenvalue))
  if (min(infodf1$total_cn[infodf1$lenvalue == max(infodf1$lenvalue)]) > 1){
    cnvdf$cnbase[cnvdf$arm == i] <- min(infodf1$total_cn[infodf1$lenvalue == max(infodf1$lenvalue)])
  }
}

cnvdf$cnstate <- "no"
for (i in 1:nrow(cnvdf)){
  if (cnvdf$copy[i] > 3*cnvdf$cnbase[i]){ # cnvdf$copy[i] > cnvdf$cnbase[i]+6 -- Size filter (cutting out <1000 sized amplicons)
    cnvdf$cnstate[i] <- "amp"
  } else if (round(cnvdf$copy[i]) < cnvdf$cnbase[i]){
    cnvdf$cnstate[i] <- "below_baseline"
  } else if (round(cnvdf$copy[i]) == cnvdf$cnbase[i]){
    cnvdf$cnstate[i] <- "baseline"
  } else if (round(cnvdf$copy[i]) == cnvdf$cnbase[i]+1){
    cnvdf$cnstate[i] <- "baseline+1"
  }
}

for (j in colnames(cnvdf)[startsWith(colnames(cnvdf), "clone_")]){
  for (i in 1:nrow(cnvdf)){
    if (cnvdf[i,j] > 3*cnvdf$cnbase[i] & cnvdf$cnstate[i] != "amp"){ # cnvdf[i,j] > cnvdf$cnbase[i]+6 -- Size filter (cutting out <1000 sized amplicons)
      cnvdf$cnstate[i] <- ifelse(!startsWith(cnvdf$cnstate[i], "amp_clone_"), paste0("amp_", j), paste0(cnvdf$cnstate[i], gsub("clone_", "_", j)))
    }
  }
}

# Annotation of amplified contiguous genomic regions
cnvdf$cgr <- "no"
cnvdf$cgr[startsWith(cnvdf$cnstate, "amp")] <- "amp"
cnvdf$acgr <- cnvdf$cgr
for (i in 1:nrow(cnvdf)){
  if (i == min(which(cnvdf$chr == cnvdf$chr[i]))){
    if (cnvdf$cgr[i] == "amp" & cnvdf$cgr[i+1] == "amp"){
      cnvdf$acgr[i] <- "amplicon_start"
    } else if (cnvdf$cgr[i] == "amp" & cnvdf$cgr[i+1] == "no"){
      cnvdf$acgr[i] <- "amplicon_itself"
    }
  } else if (i == max(which(cnvdf$chr == cnvdf$chr[i]))){
    if (cnvdf$cgr[i] == "amp" & cnvdf$cgr[i-1] == "amp"){
      cnvdf$acgr[i] <- "amplicon_end"
    } else if (cnvdf$cgr[i] == "amp" & cnvdf$cgr[i-1] == "no"){
      cnvdf$acgr[i] <- "amplicon_itself"
    }
  } else {
    if (cnvdf$cgr[i] == "amp" & cnvdf$cgr[i-1] == "no" & cnvdf$cgr[i+1] == "no"){
      cnvdf$acgr[i] <- "amplicon_itself"
    } else if (cnvdf$cgr[i] == "amp" & cnvdf$ssv[i] != "telomere" & cnvdf$cgr[i-1] == "no" & cnvdf$cgr[i+1] == "amp" | cnvdf$cgr[i] == "amp" & cnvdf$ssv[i] == "telomere" & cnvdf$cgr[i+1] == "amp"){
      cnvdf$acgr[i] <- "amplicon_start"
    } else if (cnvdf$cgr[i] == "amp" & cnvdf$esv[i] != "telomere" & cnvdf$cgr[i-1] == "amp" & cnvdf$cgr[i+1] == "no" | cnvdf$cgr[i] == "amp" & cnvdf$esv[i] == "telomere" & cnvdf$cgr[i-1] == "amp" ){
      cnvdf$acgr[i] <- "amplicon_end"
    }
  }
}

for (i in 2:(nrow(cnvdf)-1)){
  if (cnvdf$acgr[i] == "amplicon_start" & cnvdf$ssv_type[i] == "h2hINV" & nrow(isv[isv$V2 == cnvdf$chr[i] & isv$V6 <= cnvdf$start[i] & isv$V7 >= cnvdf$start[i],]) != 0){
    if ("h2hINV" %in% isv$type[isv$V2 == cnvdf$chr[i] & isv$V6 <= cnvdf$start[i] & isv$V7 >= cnvdf$start[i]]){
      k <- which(isv$V2 == cnvdf$chr[i] & isv$V6 <= cnvdf$start[i] & isv$V7 >= cnvdf$start[i] & isv$type == "h2hINV")
      if (i != 1 & i != nrow(cnvdf) & cnvdf$chr[i-1] == cnvdf$chr[i] & cnvdf$start[i-1] >= isv$V3[k] & cnvdf$start[i-1] <= isv$V4[k] & isv$V6[k] - isv$V3[k] <=5000 & cnvdf$copy[i-1] > 2*cnvdf$cnbase[i]){
        cnvdf$acgr[i] <- "amp"
        cnvdf$acgr[i-1] <- "amplicon_start"
        if (cnvdf$acgr[i-2] == "amp"){
          cnvdf$acgr[i-2] <- "amplicon_end"
          print(paste0("Amplicon extension by FBI led to conjoining amplifications in ", w, " i = ", i))
        }
      }
    }
  } else if (cnvdf$acgr[i] == "amplicon_end" & cnvdf$esv_type[i] == "t2tINV" & nrow(isv[isv$V2 == cnvdf$chr[i] & isv$V3 <= cnvdf$end[i] & isv$V4 >= cnvdf$end[i],]) != 0){
    if ("t2tINV" %in% isv$type[isv$V2 == cnvdf$chr[i] & isv$V3 <= cnvdf$end[i] & isv$V4 >= cnvdf$end[i]]){
      k <- which(isv$V2 == cnvdf$chr[i] & isv$V3 <= cnvdf$end[i] & isv$V4 >= cnvdf$end[i] & isv$type == "t2tINV")
      if (i != 1 & i != nrow(cnvdf) & cnvdf$chr[i] == cnvdf$chr[i+1] & cnvdf$end[i+1] >= isv$V6[k] & cnvdf$end[i+1] <= isv$V7[k] & isv$V6[k] - isv$V3[k] <=5000 & cnvdf$copy[i+1] > 2*cnvdf$cnbase[i]){
        cnvdf$acgr[i] <- "amp"
        cnvdf$acgr[i+1] <- "amplicon_end"
        if (cnvdf$acgr[i+2] == "amp"){
          cnvdf$acgr[i+2] <- "amplicon_start"
          print(paste0("Amplicon extension by FBI led to conjoining amplifications in ", w, " i = ", i))
        }
      }
    }
  }
}
if (sum(cnvdf$cnbase >= 6) != 0){
  print(paste0("high CN baseline (>=6): ", paste(unique(cnvdf$chr[cnvdf$cnbase >= 6]), collapse = ";")))
}

cnvdf$oncogene <- ""
for (i in which(startsWith(cnvdf$acgr, "amp"))){
  if (nrow(drvamp[drvamp$chromosome == cnvdf$chr[i] & drvamp$start <= cnvdf$end[i] & drvamp$end >= cnvdf$start[i],]) != 0){
    cnvdf$oncogene[i] <- paste0(drvamp$gene[drvamp$chromosome == cnvdf$chr[i] & drvamp$start <= cnvdf$end[i] & drvamp$end >= cnvdf$start[i]], collapse = ";")
  }
}
write.csv(cnvdf, paste0("./analysis/normalization/acgr_subclonal_cn_v4/", w, "_clusters_0.7_1000_subclonal_cn_acgr.csv"), row.names=F, quote=F)


df <- cnvdf[cnvdf$acgr %in% c("amplicon_start", "amplicon_end", "amplicon_itself"),]
if (nrow(df) == 0){
  stop("No amplification found — stopping execution.")
}

infodf <- data.frame(matrix(ncol=20, nrow=0))
colnames(infodf) <- c("study_id", "chromosome", "start", "end", "mcj_s", "mcj_e", "class_s", "svclass_s", "svname_s", "svf_s", "svlen_s", "leftcn_s", "rightcn_s", "class_e", "svclass_e", "svname_e", "svf_e", "svlen_e", "leftcn_e", "rightcn_e")
colnames(isv) <- c("id", "chr1", "start1", "end1", "chr2", "start2", "end2", "svname", "quality", "strand1", "strand2", "type")

isv$pos1 <- isv$start1
isv$pos2 <- isv$start2
isv$svlen <- NA
isv$svlen[isv$chr1 == isv$chr2] <- isv$pos2[isv$chr1 == isv$chr2] - isv$pos1[isv$chr1 == isv$chr2]

print(paste0("Processing ", w))
for (i in 1:nrow(df)){
  if (df$acgr[i] == "amplicon_itself"){
    infodf[nrow(infodf)+1,] <- NA
    infodf$study_id[nrow(infodf)] <- w
    infodf$chromosome[nrow(infodf)] <- df$chr[i]
    infodf$start[nrow(infodf)] <- df$start[i]
    infodf$end[nrow(infodf)] <- df$end[i]

    infodf$mcj_s[nrow(infodf)] <- ifelse(df$ssv[i] != "telomere", ifelse(cnvdf$cnstate[which(cnvdf$chr==df$chr[i] & cnvdf$start==df$start[i])-1] %in% c("baseline", "baseline+1", "below_baseline"), cnvdf$cnstate[which(cnvdf$chr==df$chr[i] & cnvdf$start==df$start[i])-1], "no"), "no")
    infodf$svclass_s[nrow(infodf)] <- df$ssv_type[i]
    infodf$svname_s[nrow(infodf)] <- ifelse(df$ssv[i] %in% c("telomere", "centromere", ""), "no", df$ssv[i])
    infodf$svf_s[nrow(infodf)] <- df$svf[i]
    infodf$svlen_s[nrow(infodf)] <- ifelse(df$ssv[i] %in% c("telomere", "centromere", ""), NA, isv$svlen[isv$svname == paste0(str_sub(df$ssv[i], 1, -2), "o") | isv$svname == paste0(str_sub(df$ssv[i], 1, -2), "h")])
    infodf$leftcn_s[nrow(infodf)] <- ifelse(df$ssv[i] != "telomere", cnvdf$copy[which(cnvdf$chr==df$chr[i] & cnvdf$start==df$start[i])-1], NA)
    infodf$rightcn_s[nrow(infodf)] <- df$copy[i]

    infodf$mcj_e[nrow(infodf)] <- ifelse(df$esv[i] != "telomere", ifelse(cnvdf$cnstate[which(cnvdf$chr==df$chr[i] & cnvdf$start==df$start[i])+1] %in% c("baseline", "baseline+1", "below_baseline"), cnvdf$cnstate[which(cnvdf$chr==df$chr[i] & cnvdf$start==df$start[i])+1], "no"), "no")
    infodf$svclass_e[nrow(infodf)] <- df$esv_type[i]
    infodf$svname_e[nrow(infodf)] <- ifelse(df$esv[i] %in% c("telomere", "centromere", ""), "no", df$esv[i])
    infodf$svf_e[nrow(infodf)] <- df$evf[i]
    infodf$svlen_e[nrow(infodf)] <- ifelse(df$esv[i] %in% c("telomere", "centromere", ""), NA, isv$svlen[isv$svname == paste0(str_sub(df$esv[i], 1, -2), "o") | isv$svname == paste0(str_sub(df$esv[i], 1, -2), "h")])
    infodf$leftcn_e[nrow(infodf)] <- df$copy[i]
    infodf$rightcn_e[nrow(infodf)] <- ifelse(df$esv[i] != "telomere", cnvdf$copy[which(cnvdf$chr==df$chr[i] & cnvdf$start==df$start[i])+1], NA)

  } else if (df$acgr[i] == "amplicon_start"){
    if (df$acgr[i+1] != "amplicon_end"){
      print(paste0(w, " unexpected error in acgr"))
      break
    }
    infodf[nrow(infodf)+1,] <- NA
    infodf$study_id[nrow(infodf)] <- w
    infodf$chromosome[nrow(infodf)] <- df$chr[i]
    infodf$start[nrow(infodf)] <- df$start[i]
    infodf$end[nrow(infodf)] <- df$end[i+1]

    infodf$mcj_s[nrow(infodf)] <- ifelse(df$ssv[i] != "telomere", ifelse(cnvdf$cnstate[which(cnvdf$chr==df$chr[i] & cnvdf$start==df$start[i])-1] %in% c("baseline", "baseline+1", "below_baseline"), cnvdf$cnstate[which(cnvdf$chr==df$chr[i] & cnvdf$start==df$start[i])-1], "no"), "no")
    infodf$svclass_s[nrow(infodf)] <- df$ssv_type[i]
    infodf$svname_s[nrow(infodf)] <- ifelse(df$ssv[i] %in% c("telomere", "centromere", ""), "no", df$ssv[i])
    infodf$svf_s[nrow(infodf)] <- df$svf[i]
    infodf$svlen_s[nrow(infodf)] <- ifelse(df$ssv[i] %in% c("telomere", "centromere", ""), NA, isv$svlen[isv$svname == paste0(str_sub(df$ssv[i], 1, -2), "o") | isv$svname == paste0(str_sub(df$ssv[i], 1, -2), "h")])
    infodf$leftcn_s[nrow(infodf)] <- ifelse(df$ssv[i] != "telomere", cnvdf$copy[which(cnvdf$chr==df$chr[i] & cnvdf$start==df$start[i])-1], NA)
    infodf$rightcn_s[nrow(infodf)] <- df$copy[i]

    infodf$mcj_e[nrow(infodf)] <- ifelse(df$esv[i+1] != "telomere", ifelse(cnvdf$cnstate[which(cnvdf$chr==df$chr[i+1] & cnvdf$start==df$start[i+1])+1] %in% c("baseline", "baseline+1", "below_baseline"), cnvdf$cnstate[which(cnvdf$chr==df$chr[i+1] & cnvdf$start==df$start[i+1])+1], "no"), "no")
    infodf$svclass_e[nrow(infodf)] <- df$esv_type[i+1]
    infodf$svname_e[nrow(infodf)] <- ifelse(df$esv[i+1] %in% c("telomere", "centromere", ""), "no", df$esv[i+1])
    infodf$svf_e[nrow(infodf)] <- df$evf[i+1]
    infodf$svlen_e[nrow(infodf)] <- ifelse(df$esv[i+1] %in% c("telomere", "centromere", ""), NA, isv$svlen[isv$svname == paste0(str_sub(df$esv[i+1], 1, -2), "o") | isv$svname == paste0(str_sub(df$esv[i+1], 1, -2), "h")])
    infodf$leftcn_e[nrow(infodf)] <- df$copy[i+1]
    infodf$rightcn_e[nrow(infodf)] <- ifelse(df$esv[i+1] != "telomere", cnvdf$copy[which(cnvdf$chr==df$chr[i+1] & cnvdf$start==df$start[i+1])+1], NA)

  } else if (df$acgr[i] == "amplicon_end"){

  }
}

infodf$regionlength <- infodf$end - infodf$start # will filter based on region length 50K
infodf$oncogene <- ""
infodf$peakgene <- ""
infodf$peakregion <- ""
infodf$peaklength <- NA

for (i in 1:nrow(infodf)){
  if (nrow(drvamp[drvamp$chromosome == infodf$chromosome[i] & drvamp$start <= infodf$end[i] & drvamp$end >= infodf$start[i],]) != 0){
    infodf$oncogene[i] <- paste(drvamp$gene[drvamp$chromosome == infodf$chromosome[i] & drvamp$start <= infodf$end[i] & drvamp$end >= infodf$start[i]], collapse = ";")
  }
}

k <- 0
for (i in 1:nrow(infodf)){
  k <- which(rowSums(cnvdf[,c("copy", colnames(cnvdf)[startsWith(colnames(cnvdf), "clone_")])] == max(cnvdf[cnvdf$chr == infodf$chromosome[i] & cnvdf$start >= infodf$start[i] & cnvdf$end <= infodf$end[i],c("copy", colnames(cnvdf)[startsWith(colnames(cnvdf), "clone_")])])) > 0)[1]
  infodf$peakregion[i] <- paste0(cnvdf$chr[k], "_", cnvdf$start[k], "_", cnvdf$end[k])
  infodf$peaklength[i] <- cnvdf$end[k] - cnvdf$start[k]
  if (length(genedf$gene[genedf$chromosome == cnvdf$chr[k] & genedf$start <= cnvdf$end[k] & genedf$end >= cnvdf$start[k]][genedf$gene[genedf$chromosome == cnvdf$chr[k] & genedf$start <= cnvdf$end[k] & genedf$end >= cnvdf$start[k]] %in% drvamp$gene]) != 0){
    infodf$peakgene[i] <- genedf$gene[genedf$chromosome == cnvdf$chr[k] & genedf$start <= cnvdf$end[k] & genedf$end >= cnvdf$start[k]][genedf$gene[genedf$chromosome == cnvdf$chr[k] & genedf$start <= cnvdf$end[k] & genedf$end >= cnvdf$start[k]] %in% drvamp$gene][1]
  } else {
    infodf$peakgene[i] <- genedf$gene[genedf$chromosome == cnvdf$chr[k] & genedf$start <= cnvdf$end[k] & genedf$end >= cnvdf$start[k]][1]
  }
}

infodf$measureregion <- ""
for (i in which(infodf$regionlength > 50000 & infodf$oncogene != "")){
  for (j in strsplit(infodf$oncogene[i], ";", fixed=T)[[1]]){
    tempdf <- cnvdf[cnvdf$chr == drvamp$chromosome[drvamp$gene == j] & cnvdf$start <= drvamp$end[drvamp$gene == j] & cnvdf$end >= drvamp$start[drvamp$gene == j],]
    if (nrow(tempdf) != 0 && max(tempdf$end) - min(tempdf$start) + 1 > 50000){
      chri <- drvamp$chromosome[drvamp$gene == j]
      starti <- min(tempdf$start)
      endi <- max(tempdf$end)
    } else {
      tempdf <- cnvdf[cnvdf$chr == infodf$chromosome[i] & cnvdf$start <= infodf$end[i] & cnvdf$end >= infodf$start[i],]
      chri <- drvamp$chromosome[drvamp$gene == j]
      if (drvamp$start[drvamp$gene == j] - 25000 > min(tempdf$start) & drvamp$end[drvamp$gene == j] + 25000 < max(tempdf$end)){
        starti <- max(tempdf$start[tempdf$start <= drvamp$start[drvamp$gene == j] - 25000])
        endi <- min(tempdf$end[tempdf$end >= drvamp$end[drvamp$gene == j] + 25000])
      } else if (drvamp$start[drvamp$gene == j] - 50000 > min(tempdf$start) & drvamp$end[drvamp$gene == j] < max(tempdf$end)){
        starti <- max(tempdf$start[tempdf$start <= drvamp$start[drvamp$gene == j] - 50000])
        endi <- min(tempdf$end[tempdf$end >= drvamp$end[drvamp$gene == j]])
      } else if (drvamp$start[drvamp$gene == j] > min(tempdf$start) & drvamp$end[drvamp$gene == j] + 50000 < max(tempdf$end)){
        starti <- max(tempdf$start[tempdf$start <= drvamp$start[drvamp$gene == j]])
        endi <- min(tempdf$end[tempdf$end >= drvamp$end[drvamp$gene == j] + 50000])
      } else {
        chri <- drvamp$chromosome[drvamp$gene == j]
        starti <- drvamp$start[drvamp$gene == j]
        endi <- drvamp$end[drvamp$gene == j]
        print(paste0("Amplicon does not fit for analysis for ", j))
      }
      
    }
    infodf$measureregion[i] <- ifelse(infodf$measureregion[i] == "", paste0(j, "_", chri, "_", starti, "_", endi), paste0(infodf$measureregion[i], paste0(";", j, "_", chri, "_", starti, "_", endi)))
  }
}
infodf$regionid <- ""
infodf$regionid[infodf$measureregion != ""] <- paste0(w, "_region", c(1:length(infodf$regionid[infodf$measureregion != ""])))
write.csv(infodf, paste0("./analysis/normalization/amplified_region_v4/", w, ".amplified_region.csv"), row.names = F, quote = F)

# Annotate per-cell CN for all relevant loci
df <- fread(paste0("./analysis/sv-seg-cn-swap/results/", w, "/copy_number/", w, ".normalized_merged.csv.gz"))
df[, length := end - start]
qc <- fread(paste0("./analysis/metadata/qc/", w, "/", w, "_qc.csv.gz"))
cndf <- qc[
  quality >= 0.75 &
    total_mapped_reads > 250000 &
    is_tumor_cell == "yes" &
    is_control == F &
    is_contaminated == F &
    is_s_phase == F &
    !grepl("CONTROL", sample_type),
  .(cell_id, sample_id, sample_type, quality, total_mapped_reads)
]

regions <- unique(strsplit(paste(infodf$measureregion[infodf$measureregion != ""], collapse = ";"), ";", fixed = TRUE)[[1]])
if (length(regions) != 0){
  regiondt <- rbindlist(lapply(regions, function(x) {
    parts <- strsplit(x, "_", fixed = TRUE)[[1]]
    data.table(
      gene = parts[1],
      chr = parts[2],
      start = as.numeric(parts[3]),
      end = as.numeric(parts[4])
    )
  }))
  
  setkey(df, chr, start, end)
  setkey(regiondt, chr, start, end)
  
  # Perform overlap join
  ol <- foverlaps(df, regiondt, by.x = c("chr", "start", "end"), by.y = c("chr", "start", "end"), type = "within", nomatch = 0)
  oncogene_cn <- ol[
    , .(weighted_copy = weighted.mean(copy, length, na.rm = TRUE)),
    by = .(cell_id, gene)
  ]
  
  oncogene_cn_wide <- dcast(oncogene_cn, cell_id ~ gene, value.var = "weighted_copy")
  
  cndf <- merge(cndf, oncogene_cn_wide, by = "cell_id", all.x = TRUE)
  fwrite(cndf, paste0("./analysis/normalization/oncogene_cn_v4/", w, ".sv-seg-cn.oncogene.cn.csv"))
}


