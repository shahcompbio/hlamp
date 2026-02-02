# Part 1: filtering SVs with no read support
ptlist <- indi$individual[indi$normal_wgs == "yes"]
ptlist <- ptlist[!(ptlist %in% c("SA535", "SA609"))]
infoline <- c()
for (i in 1:length(ptlist)){
  infoline[i] <- paste0('zcat ./analysis/hmf/gridss_somatic/', ptlist[i], '/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vcf.gz | grep -v "VF=0" > ./analysis/hmf/gridss_somatic/', ptlist[i], '/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.vcf &&')
}
writeLines(c('#!/bin/bash', infoline), "./analysis/hmf/gridss_somatic/additional_filter.sh")

# Part 2: VCF to BEDPE
library(StructuralVariantAnnotation)
library(rtracklayer)

for (i in which(indi$individual != "2765_2")){
  w <- indi$individual[i]
  print(paste0("Processing ", w))

  vcf = readVcf(paste0("analysis/hmf/gridss_", ifelse(indi$normal_wgs[indi$individual == w] == "yes", "somatic", "tumoronly"), "/", w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.vcf"))

  # Export breakpoints to BEDPE
  vcf <- vcf[ fixed(vcf)$FILTER == "PASS" ]
  bpgr <- breakpointRanges(vcf)
  write.table(breakpointgr2bedpe(bpgr), file=paste0("analysis/hmf/gridss_", ifelse(indi$normal_wgs[indi$individual == w] == "yes", "somatic", "tumoronly"), "/", w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.paired.bedpe"), sep="\t", quote=FALSE, col.names=FALSE)

  # Export single breakends to BED
  begr = breakendRanges(vcf)
  begr$score = begr$QUAL
  export(begr, con=paste0("analysis/hmf/gridss_", ifelse(indi$normal_wgs[indi$individual == w] == "yes", "somatic", "tumoronly"), "/", w, "/COMPOSITE_TUMOR_LIBRARIES.gripss.filtered.vf.singlebreakend.bed"))
}

