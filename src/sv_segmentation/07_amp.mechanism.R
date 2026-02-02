library(stringr)

# Analyzing amplified regions and their mechanisms
ptlist <- gsub(".amplified_region.csv", "", list.files(path = "./analysis/normalization/amplified_region_v4", pattern = ".amplified_region.csv$", recursive = F, full.names = F))
for (i in 1:length(ptlist)){
  if (i == 1){
    df <- read.csv(paste0("./analysis/normalization/amplified_region_v4/", ptlist[i], ".amplified_region.csv"))
  } else {
    tempdf <- read.csv(paste0("./analysis/normalization/amplified_region_v4/", ptlist[i], ".amplified_region.csv"))
    df <- rbind(df, tempdf)
  }
}
colnames(df)[1] <- "individual"

# Annotation of SV clusters
df$svcluster_s <- ""
df$svcluster_e <- ""
for (w in ptlist){
  print(paste0("Processing ", w))
  isv <- read.csv(paste0("./analysis/clustering/clusters_gridss_090325/", w, ".sv_clusters_and_footprints.tsv"), header=F, as.is=T, sep="\t")
  isv <- isv[isv$V12 != 1,]

  for (i in which(df$individual == w)){
    chri <- df$chromosome[i]
    starti <- df$start[i]
    endi <- df$end[i]

    df$svcluster_s[i] <- ifelse(nrow(isv[str_sub(isv$V7, 1, -2) == str_sub(df$svname_s[i], 1, -2),]) != 0, paste0("cluster", isv$V11[str_sub(isv$V7, 1, -2) == str_sub(df$svname_s[i], 1, -2)], "[", isv$V12[str_sub(isv$V7, 1, -2) == str_sub(df$svname_s[i], 1, -2)], "]"), "")
    df$svcluster_e[i] <- ifelse(nrow(isv[str_sub(isv$V7, 1, -2) == str_sub(df$svname_e[i], 1, -2),]) != 0, paste0("cluster", isv$V11[str_sub(isv$V7, 1, -2) == str_sub(df$svname_e[i], 1, -2)], "[", isv$V12[str_sub(isv$V7, 1, -2) == str_sub(df$svname_e[i], 1, -2)], "]"), "")
  }
}

# Chromothripsis annotation
df$chromothripsis <- ""
for (i in which(df$svcluster_s != "")){
  print(i)
  w <- df$individual[i]
  isv <- read.csv(paste0("./analysis/clustering/clusters_gridss_090325/", w, ".sv_clusters_and_footprints.tsv"), header=F, as.is=T, sep="\t")
  isv <- isv[isv$V11 == as.numeric(strsplit(gsub("cluster", "", df$svcluster_s[i]), "[", fixed=T)[[1]][1]),]
  isv <- isv[!(isv$V1 == isv$V4 & isv$V5 - isv$V2 < 100 & isv$V9 == "+" & isv$V10 == "-"),]
  isv$type <- paste0(isv$V9, isv$V10)
  for (j in names(sort(table(c(isv$V1[isv$V1 == isv$V4], isv$V4[isv$V1 == isv$V4]))[table(c(isv$V1[isv$V1 == isv$V4], isv$V4[isv$V1 == isv$V4])) >= 20], decreasing = T))){
    observed_counts <- unname(table(isv$type[isv$V1 == j & isv$V4 == j]))
    if (length(observed_counts) == 4 && permutation_test(observed_counts)$p_value >= 0.05 && nrow(isv[isv$V1 == j & isv$V4 == j,]) >= 10){
      #acdf <- read.csv(paste0("analysis/segs/", ifelse(indi$normal_wgs[indi$individual == w] == "yes", "somatic", "tumoronly"), "/", w, "/cnv.acgr.tsv"), header=T, as.is=T, sep="\t")
      acdf <- read.csv(paste0("./analysis/normalization/acgr_subclonal_cn_v4/", w, "_clusters_0.7_1000_subclonal_cn_acgr.csv"))
      acdf <- acdf[acdf$chr == j & acdf$start >= min(c(isv$V2[isv$V1 == j], isv$V5[isv$V4 == j])) & acdf$end <= max(c(isv$V6[isv$V4 == j], isv$V3[isv$V1 == j])),]
      if (sum(grepl("amp", acdf$acgr)[-1] != grepl("amp", acdf$acgr)[-length(grepl("amp", acdf$acgr))]) >= 10){
        #print(paste0(w, " has chromothripsis on the s cluster, in chromosome ", j))
        if (df$chromothripsis[i] == ""){
          df$chromothripsis[i] <- paste0("chr", j, "[", nrow(isv[isv$V1 == j & isv$V4 == j,]), "]")
        } else {
          df$chromothripsis[i] <- paste(unique(c(strsplit(df$chromothripsis[i], ";", fixed = T)[[1]], paste0("chr", j, "[", nrow(isv[isv$V1 == j & isv$V4 == j,]), "]"))), collapse = ";")
        }
      }
    }
  }
}
for (i in which(df$svcluster_e != "")){
  print(i)
  w <- df$individual[i]
  isv <- read.csv(paste0("./analysis/clustering/clusters_gridss_090325/", w, ".sv_clusters_and_footprints.tsv"), header=F, as.is=T, sep="\t")
  isv <- isv[isv$V11 == as.numeric(strsplit(gsub("cluster", "", df$svcluster_e[i]), "[", fixed=T)[[1]][1]),]
  isv <- isv[!(isv$V1 == isv$V4 & isv$V5 - isv$V2 < 100 & isv$V9 == "+" & isv$V10 == "-"),]
  isv$type <- paste0(isv$V9, isv$V10)
  for (j in names(sort(table(c(isv$V1[isv$V1 == isv$V4], isv$V4[isv$V1 == isv$V4]))[table(c(isv$V1[isv$V1 == isv$V4], isv$V4[isv$V1 == isv$V4])) >= 20], decreasing = T))){
    observed_counts <- unname(table(isv$type[isv$V1 == j & isv$V4 == j]))
    if (length(observed_counts) == 4 && permutation_test(observed_counts)$p_value >= 0.05 && nrow(isv[isv$V1 == j & isv$V4 == j,]) >= 10){
      #acdf <- read.csv(paste0("analysis/segs/", ifelse(indi$normal_wgs[indi$individual == w] == "yes", "somatic", "tumoronly"), "/", w, "/cnv.acgr.tsv"), header=T, as.is=T, sep="\t")
      acdf <- read.csv(paste0("./analysis/normalization/acgr_subclonal_cn_v4/", w, "_clusters_0.7_1000_subclonal_cn_acgr.csv"))
      acdf <- acdf[acdf$chr == j & acdf$start >= min(c(isv$V2[isv$V1 == j], isv$V5[isv$V4 == j])) & acdf$end <= max(c(isv$V6[isv$V4 == j], isv$V3[isv$V1 == j])),]
      if (sum(grepl("amp", acdf$acgr)[-1] != grepl("amp", acdf$acgr)[-length(grepl("amp", acdf$acgr))]) >= 10){
        #print(paste0(w, " has chromothripsis on the s cluster, in chromosome ", j))
        if (df$chromothripsis[i] == ""){
          df$chromothripsis[i] <- paste0("chr", j, "[", nrow(isv[isv$V1 == j & isv$V4 == j,]), "]")
        } else {
          df$chromothripsis[i] <- paste(unique(c(strsplit(df$chromothripsis[i], ";", fixed = T)[[1]], paste0("chr", j, "[", nrow(isv[isv$V1 == j & isv$V4 == j,]), "]"))), collapse = ";")
        }
      }
    }
  }
}

df$mech <- ""
df$mech[df$chromothripsis != ""] <- "CTX"
df$mech[(df$svclass_s == "h2hINV" & df$svlen_s <= 5000) | (df$svclass_e == "t2tINV" & df$svlen_e <= 5000)] <- "FBI"
df$mech[df$svclass_s == "DUP" & df$svclass_e == "DUP" & str_sub(df$svname_s, 1, -2) == str_sub(df$svname_e, 1, -2)] <- "TD"
df$mech[(df$svclass_s == "SBE" & df$svclass_e == "SBE") | (df$svclass_s == "SBE" & df$svclass_e == "") | (df$svclass_s == "" & df$svclass_e == "SBE") | (df$svclass_s == "" & df$svclass_e == "")] <- "UNK"
df$mech[df$individual == "2765_2"] <- "TD"
df$mech[df$mech == "" & (df$svclass_s == "TRA" | df$svclass_e == "TRA")] <- "TRA"
df$mech[df$mech == "TRA" & ((df$svclass_s == "TRA" & df$svcluster_s != "") | (df$svclass_e == "TRA" & df$svcluster_e != ""))] <- "CGR-TRA"
df$mech[df$mech == ""] <- "CGR-NOS"


for (i in which(df$mech %in% c("CGR-NOS", "UNK"))){
  w <- df$individual[i]

  isv <- read.csv(paste0("analysis/hmf/gridss_", ifelse(w %in% indi$individual[indi$normal_wgs == "yes"], "somatic", "tumoronly"), "/", w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.annotated.bedpe"), header=T, as.is=T, sep="\t")
  if (w %in% c("SPECTRUM-OV-004", "SPECTRUM-OV-022", "SPECTRUM-OV-046", "SPECTRUM-OV-129")){
    isv <- isv[!(isv$type == "TRA" & isv$vf <= 4),]
  }
  isv$type[isv$type %in% c("h2hFBI", "t2tFBI")] <- "FBI"

  chri <- df$chromosome[i]
  posi <- round((df$start[i] + df$end[i])/2)
  isv <- isv[(isv$chr1 == chri & isv$start1 >= posi - 5000000 & isv$end1 <= posi + 5000000) | (isv$chr2 == chri & isv$start2 >= posi - 5000000 & isv$end2 <= posi + 5000000),]
  if (nrow(isv) >= 10 && names(sort(table(isv$type), decreasing = T))[1] == "FBI"){
    df$mech[i] <- "CGR-FBI"

  }
}
ampdf <- df

options(scipen = 999)
ampdf$mean_cn <- NA
ampdf$sd_cn <- NA
ampdf$cn10 <- NA

for (w in unique(ampdf$individual)){
  print(paste0("Processing ", w))
  cnvdf <- read.csv(paste0("./analysis/normalization/region_cn_final/", w, ".sv-seg-cn.region.cn.csv"))
  for (i in which(ampdf$individual == w)){
    genei <- paste0("chr", ampdf$chromosome[i], "_", ampdf$start[i], "_", ampdf$end[i])
    ampdf$mean_cn[i] <- mean(cnvdf[,genei], na.rm = T)
    ampdf$sd_cn[i] <- sd(cnvdf[,genei], na.rm = T)
    ampdf$cn10[i] <- sum(cnvdf[,genei] >= 10)/nrow(cnvdf)
    rm(genei)
  }
}
write.table(ampdf, "amplification.info.tsv", row.names=F, col.names=T, quote=F, sep="\t")

permutation.test <- function(observed_counts, num_permutations = 1000, seed = 123) {
  if (length(observed_counts) != 4) {
    stop("Observed counts should be a vector of four numbers.")
  }
  
  total_count <- sum(observed_counts)
  
  permutation_counts <- matrix(0, nrow = num_permutations, ncol = 4)
  
  set.seed(seed)
  
  for (i in 1:num_permutations) {
    random_counts <- sample(c("++", "+-", "-+", "--"), total_count, replace = TRUE)
    permutation_counts[i, ] <- table(factor(random_counts, levels = c("++", "+-", "-+", "--")))
  }
  
  expected_count <- total_count / 4
  observed_statistic <- sum((observed_counts - expected_count)^2 / expected_count)
  
  permutation_statistics <- apply(permutation_counts, 1, function(x) {
    sum((x - expected_count)^2 / expected_count)
  })
  
  p_value <- mean(permutation_statistics >= observed_statistic)
  
  return(list(observed_statistic = observed_statistic, p_value = p_value))
}




