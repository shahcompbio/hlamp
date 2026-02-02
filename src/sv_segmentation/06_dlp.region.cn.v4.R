library(dplyr)
library(stringr)
library(signals)
library(data.table)
library(anndata)

options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
w <- args

# Load region information
infodf <- read.csv(paste0("./analysis/normalization/amplified_region_v4/", w, ".amplified_region.csv"))

# Annotate per-cell CN for all relevant loci
if (w %in% c("2765_2")){
  df <- fread(paste0("./analysis/sv-seg-cn-mouse/results/", w, "/copy_number/", w, ".normalized_merged.csv.gz"))
} else {
  df <- fread(paste0("./analysis/sv-seg-cn-swap/results/", w, "/copy_number/", w, ".normalized_merged.csv.gz"))
}
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

if (nrow(infodf) != 0){
  regiondt <- data.table(
    region <- paste0("chr", infodf$chromosome, "_", infodf$start, "_", infodf$end),
    chr <- infodf$chromosome,
    start <- infodf$start,
    end <- infodf$end
  )
  colnames(regiondt) <- c("region", "chr", "start", "end")
  regiondt$chr <- as.character(regiondt$chr)
  
  setkey(df, chr, start, end)
  setkey(regiondt, chr, start, end)
  
  # Perform overlap join
  ol <- foverlaps(df, regiondt, by.x = c("chr", "start", "end"), by.y = c("chr", "start", "end"), type = "within", nomatch = 0)
  region_cn <- ol[
    , .(weighted_copy = weighted.mean(copy, length, na.rm = TRUE)),
    by = .(cell_id, region)
  ]
  
  region_cn_wide <- dcast(region_cn, cell_id ~ region, value.var = "weighted_copy")
  
  cndf <- merge(cndf, region_cn_wide, by = "cell_id", all.x = TRUE)
  fwrite(cndf, paste0("./analysis/normalization/region_cn_v4/", w, ".sv-seg-cn.region.cn.csv"))
}

