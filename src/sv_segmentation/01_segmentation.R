library(dplyr)
library(stringr)
library(signals)

ptlist <- indi$individual[indi$normal_wgs == "yes"]

binsize <- 2000000
options(scipen = 999)
for (w in ptlist){
  print(paste0("Processing ", w))

  if (!dir.exists(paste0("./analysis/segs/somatic/", w))){
    dir.create(paste0("./analysis/segs/somatic/", w))
  }

  df <- read.csv(paste0("./analysis/hmf/gridss_somatic/", w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.vcf"), header=F, as.is=T, sep="\t", comment.char = "#")
  df <- df[df$V7 == "PASS",]
  df$strand <- ""
  df$vf <- NA
  for (i in 1:nrow(df)){
    if (str_sub(df$V5[i], 1, 1) %in% c("]", "[")){
      df$strand[i] <- "-"
    } else {
      df$strand[i] <- "+"
    }
    df$vf[i] <- as.numeric(strsplit(df[i,ncol(df)-2], ":", fixed = T)[[1]][length(strsplit(df[i,ncol(df)-2], ":", fixed = T)[[1]])])
  }
  df <- df[,c("V1", "V2", "V3", "strand", "vf")]
  df <- df[df$V1 %in% c(c(1:22), "X", "Y"),]

  segs <- as.data.frame(matrix(NA, nrow = 0, ncol = 8))
  colnames(segs) <- c("chr", "start", "end", "id", "ssv", "esv", "svf", "evf")
  for (i in 1:ncol(segs)){
    if (colnames(segs)[i]  %in% c("chr", "id", "ssv", "esv")){
      segs[,i] <- as.character(segs[,i])
    } else {
      segs[,i] <- as.numeric(segs[,i])
    }
  }
  for (i in as.character(c(c(1:22), "X", "Y"))){
    tempdf <- as.data.frame(matrix(NA, nrow = 0, ncol = 8))
    colnames(tempdf) <- c("chr", "start", "end", "id", "ssv", "esv", "svf", "evf")

    for (j in 1:(nrow(df[df$V1 == i,])+1)){
      tempdf[nrow(tempdf)+1,] <- c(i, NA, NA, paste0("chr", i, "-", j), "", "", NA, NA)
    }
    tempdf$start <- as.numeric(tempdf$start)
    tempdf$end <- as.numeric(tempdf$end)
    tempdf$svf <- as.numeric(tempdf$svf)
    tempdf$evf <- as.numeric(tempdf$evf)

    if (nrow(df[df$V1 == i,]) != 0){
      for (j in 1:nrow(df[df$V1 == i,])){
        if (df$strand[df$V1 == i][j] == "+"){
          tempdf$end[j] <- df$V2[df$V1 == i][j]
          tempdf$esv[j] <- df$V3[df$V1 == i][j]
          tempdf$evf[j] <- df$vf[df$V1 == i][j]
        } else {
          tempdf$start[j+1] <- df$V2[df$V1 == i][j]
          tempdf$ssv[j+1] <- df$V3[df$V1 == i][j]
          tempdf$svf[j+1] <- df$vf[df$V1 == i][j]
        }
      }
      tempdf$start[1] <- 1
      tempdf$ssv[1] <- "telomere"
      tempdf$end[nrow(tempdf)] <- hg19$length[hg19$chr == i]
      tempdf$esv[nrow(tempdf)] <- "telomere"
      for (j in 1:nrow(tempdf)){
        if (is.na(tempdf$start[j])){
          tempdf$start[j] <- tempdf$end[j-1]+1
        }
        if (is.na(tempdf$end[j])){
          tempdf$end[j] <- tempdf$start[j+1]-1
        }
      }
      tempdf <- tempdf[tempdf$end >= tempdf$start,]
    } else {
      tempdf$start[1] <- 1
      tempdf$ssv[1] <- "telomere"
      tempdf$end[nrow(tempdf)] <- hg19$length[hg19$chr == i]
      tempdf$esv[nrow(tempdf)] <- "telomere"
    }

    # centromere start
    if (!(hg19$centro_start[hg19$chr == i] %in% tempdf$start) & !(hg19$centro_start[hg19$chr == i] %in% tempdf$end)){
      imp_num <- which(tempdf$start < hg19$centro_start[hg19$chr == i] & tempdf$end > hg19$centro_start[hg19$chr == i])
      tempdf[nrow(tempdf)+1,] <- tempdf[imp_num,]
      tempdf$end[imp_num] <- hg19$centro_start[hg19$chr == i]
      tempdf$esv[imp_num] <- ""
      tempdf$evf[imp_num] <- NA
      tempdf$start[nrow(tempdf)] <- hg19$centro_start[hg19$chr == i]+1
      tempdf$ssv[nrow(tempdf)] <- "centromere"
      tempdf$svf[nrow(tempdf)] <- NA
      tempdf <- tempdf[order(tempdf$start, decreasing = F),]
    }

    # centromere end
    if (!(hg19$centro_end[hg19$chr == i] %in% tempdf$start) & !(hg19$centro_end[hg19$chr == i] %in% tempdf$end)){
      imp_num <- which(tempdf$start < hg19$centro_end[hg19$chr == i] & tempdf$end > hg19$centro_end[hg19$chr == i])
      tempdf[nrow(tempdf)+1,] <- tempdf[imp_num,]
      tempdf$end[imp_num] <- hg19$centro_end[hg19$chr == i]
      tempdf$esv[imp_num] <- "centromere"
      tempdf$evf[imp_num] <- NA
      tempdf$start[nrow(tempdf)] <- hg19$centro_end[hg19$chr == i]+1
      tempdf$ssv[nrow(tempdf)] <- ""
      tempdf$svf[nrow(tempdf)] <- NA
      tempdf <- tempdf[order(tempdf$start, decreasing = F),]
    }

    tempdf$seglen <- tempdf$end - tempdf$start + 1

    for (j in which(tempdf$seglen >= binsize*2 & tempdf$ssv != "centromere" & tempdf$esv != "centromere")){
      for (k in c((ceiling(tempdf$end[j]/binsize)-1):ceiling(tempdf$start[j]/binsize))){
        tempdf[nrow(tempdf)+1,] <- tempdf[j,]
        tempdf$end[j] <- k*binsize
        tempdf$esv[j] <- ""
        tempdf$evf[j] <- NA
        tempdf$start[nrow(tempdf)] <- k*binsize + 1
        tempdf$ssv[nrow(tempdf)] <- ""
        tempdf$svf[nrow(tempdf)] <- NA
      }
    }
    tempdf <- tempdf[order(tempdf$start, decreasing = F),]
    if (nrow(tempdf[tempdf$start > tempdf$end,]) != 0){
      print(paste0("chromosome ", i))
    }
    tempdf <- tempdf[tempdf$start <= tempdf$end,]
    tempdf$id <- paste0("chr", i, "-", 1:nrow(tempdf))

    segs <- rbind(segs, tempdf[,colnames(tempdf)[colnames(tempdf) != "seglen"]])
  }
  segs$seglen <- segs$end - segs$start + 1
  write.table(segs, paste0("./analysis/segs/somatic/", w, "/segment.sv.tsv"), row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(segs[,c(1:4)], paste0("./analysis/segs/somatic/", w, "/segment.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
}
