# VF-annotated bedpe files
for (w in ptlist){
  print(paste0("Processing ", w))
  
  svdf <- read.csv(paste0("analysis/hmf/gridss_", ifelse(indi$normal_wgs[indi$individual == w] == "yes", "somatic", "tumoronly"), "/", w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.bedpe"), header=F, as.is=T, sep="\t")
  #svdf <- read.csv(paste0("analysis/hmf/gridss_", ifelse(w %in% c("P-0042383", "P-0056027", "P-0114502"), "somatic", "tumoronly"), "/", w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.bedpe"), header=F, as.is=T, sep="\t")
  svdf <- svdf[,colnames(svdf)[colnames(svdf) != "V1"]]
  svdf$svlen <- round((svdf$V6 + svdf$V7)/2 - (svdf$V3 + svdf$V4)/2)
  svdf$svlen[svdf$V2 != svdf$V5] <- 0
  svdf$type <- ""
  svdf$type[svdf$V10 == "+" & svdf$V11 == "-"] <- "DEL"
  svdf$type[svdf$V10 == "-" & svdf$V11 == "+"] <- "DUP"
  svdf$type[svdf$V10 == "+" & svdf$V11 == "+"] <- "t2tINV"
  svdf$type[svdf$V10 == "-" & svdf$V11 == "-"] <- "h2hINV"
  svdf$type[svdf$V2 != svdf$V5] <- "TRA"
  colnames(svdf) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "svname", "qual", "ori1", "ori2", "svlen", "type")
  svdf$svname <- str_sub(svdf$svname, 1, -2)
  svdf <- svdf[!(svdf$type == "DEL" & svdf$svlen < 100),]
  
  vcf <- readVcf(paste0("analysis/hmf/gridss_", ifelse(indi$normal_wgs[indi$individual == w] == "yes", "somatic", "tumoronly"), "/", w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.vcf"))
  vcf <- vcf[fixed(vcf)$FILTER == "PASS"]
  vcf_event <- info(vcf)$EVENT
  vcf_vf <- info(vcf)$VF
  vf_map <- setNames(vcf_vf, vcf_event)
  svdf$vf <- vf_map[svdf$svname]
  
  isv <- read.csv(paste0("./analysis/clustering/clusters_gridss_090325/", w, ".sv_clusters_and_footprints.tsv"), header=F, as.is=T, sep="\t")
  isv <- isv[isv$V12 != 1,]
  isv <- isv[!(isv$V1 == isv$V4 & isv$V5 - isv$V2 < 100 & isv$V9 == "+" & isv$V10 == "-"),]
  svdf$cluster <- ""
  svdf$color <- ""
  imp_num <- 0
  for (i in as.numeric(names(table(isv$V11)[table(isv$V11) >= 2]))){
    imp_num <- imp_num + 1
    svdf$cluster[svdf$svname %in% str_sub(isv$V7[isv$V11 == i], 1, -2)] <- paste0("cluster", i, "[", nrow(isv[isv$V11 == i,]), "]")
    svdf$color[svdf$svname %in% str_sub(isv$V7[isv$V11 == i], 1, -2)] <- paste0("color=(", as.numeric(col2rgb(viridisLite::cividis(length(names(table(isv$V11)[table(isv$V11) >= 2]))))[,imp_num][1]), ",", as.numeric(col2rgb(viridisLite::cividis(length(names(table(isv$V11)[table(isv$V11) >= 2]))))[,imp_num][2]), ",", as.numeric(col2rgb(viridisLite::cividis(length(names(table(isv$V11)[table(isv$V11) >= 2]))))[,imp_num][3]), ",0.4)")
  }
  svdf$color[svdf$color == "" & svdf$chr1 != svdf$chr2] <- "color=(85,26,139,0.4)"
  svdf$color[svdf$color == "" & svdf$chr1 == svdf$chr2 & svdf$ori1 == "+" & svdf$ori2 == "-"] <- "color=(0,0,255,0.4)"
  svdf$color[svdf$color == "" & svdf$chr1 == svdf$chr2 & svdf$ori1 == "-" & svdf$ori2 == "+"] <- "color=(0,128,128,0.4)"
  svdf$color[svdf$color == "" & svdf$chr1 == svdf$chr2 & svdf$ori1 == "+" & svdf$ori2 == "+" & svdf$svlen >= 5000] <- "color=(128,128,0,0.4)"
  svdf$color[svdf$color == "" & svdf$chr1 == svdf$chr2 & svdf$ori1 == "+" & svdf$ori2 == "+" & svdf$svlen < 5000] <- "color=(139,10,80,1)"
  svdf$color[svdf$color == "" & svdf$chr1 == svdf$chr2 & svdf$ori1 == "-" & svdf$ori2 == "-" & svdf$svlen >= 5000] <- "color=(220,20,60,0.4)"
  svdf$color[svdf$color == "" & svdf$chr1 == svdf$chr2 & svdf$ori1 == "-" & svdf$ori2 == "-" & svdf$svlen < 5000] <- "color=(139,10,80,1)"
  
  svdf$type[svdf$chr1 == svdf$chr2 & svdf$ori1 == "+" & svdf$ori2 == "+" & svdf$svlen < 5000] <- "t2tFBI"
  svdf$color[svdf$chr1 == svdf$chr2 & svdf$ori1 == "+" & svdf$ori2 == "+" & svdf$svlen < 5000] <- "color=(139,10,80,1)"
  svdf$type[svdf$chr1 == svdf$chr2 & svdf$ori1 == "-" & svdf$ori2 == "-" & svdf$svlen < 5000] <- "h2hFBI"
  svdf$color[svdf$chr1 == svdf$chr2 & svdf$ori1 == "-" & svdf$ori2 == "-" & svdf$svlen < 5000] <- "color=(139,10,80,1)"
  write.table(svdf, paste0("analysis/hmf/gridss_", ifelse(w %in% indi$individual[indi$normal_wgs == "yes"], "somatic", "tumoronly"), "/", w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.annotated.bedpe"), row.names = F, col.names = T, quote = F, sep = "\t") 
}
