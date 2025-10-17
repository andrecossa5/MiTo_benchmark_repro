# Prep redeem mouse data
library(tidyverse)
library(redeemR)


##


# Utils
make_dic_ATAC <- function() {
  data(ATACWhite)
  data(RNAWhite)
  dic <- ATACWhite$V1
  names(dic) <- as.character(RNAWhite$V1)
  dic
}


##


to_ATAC <- function(name, Dic_ATAC = NULL, bclength = 16, from = c(1, 2, 3), to = c(1, 2, 3), PostFix = F) {
  
  stopifnot(is.character(name), length(from) == length(to))
  L <- nchar(name)
  NakeName <- substr(name, 1L, bclength)
  post <- ifelse(L >= bclength + 2L, substr(name, bclength + 2L, L), "")
  
  if (length(from)) {
    post_map <- setNames(to, from)
    mi <- match(post, names(post_map))
    hit <- !is.na(mi)
    post[hit] <- unname(post_map[mi[hit]])
  }
  
  new_prefix <- unname(Dic_ATAC[NakeName])
  
  if (isTRUE(PostFix)) {
    ifelse(nzchar(post), paste(new_prefix, post, sep = "_"), new_prefix)
  } else {
    new_prefix
  }
}


##


make_dic_RNA <- function() {
  data(ATACWhite)
  data(RNAWhite)
  dic <- RNAWhite$V1
  names(dic) <- as.character(ATACWhite$V1)
  dic
}


##


to_RNA <- function(name, Dic_RNA = NULL, bclength = 16,
                   from = c(1, 2, 3), to = c(1, 2, 3),
                   PostFix = TRUE) {
  
  stopifnot(is.character(name), length(from) == length(to))
  
  L <- nchar(name)
  NakeName <- substr(name, 1L, bclength)
  post <- ifelse(L >= bclength + 2L, substr(name, bclength + 2L, L), "")
  
  if (length(from)) {
    post_map <- setNames(to, from)
    mi <- match(post, names(post_map))
    hit <- !is.na(mi)
    post[hit] <- unname(post_map[mi[hit]])
  }
  
  new_prefix <- unname(Dic_RNA[NakeName])
  
  if (isTRUE(PostFix)) {
    ifelse(nzchar(post), paste(new_prefix, post, sep = "-"), new_prefix)
  } else {
    new_prefix
  }
}


##


# Read data

# Path
path <- '/Users/IEO5505/Desktop/MI_TO/last_data/redeem_mouse/GSE259284_RAW'

# Substitute below code (redeemR.read.trim()) to handle coverage data in .txt format (lists, as RDS file)

# Create trimmed-ends basecall table
# VariantsGTSummary <- redeemR.read.trim(
#  path_input, thr="T", Processed=F, rdsname="/VariantsGTSummary.S.trim5_binom.rds", edge_trim=5
# )

# Set params
edge_trim = 4
thr = "T"

# Read provided files (GEO: GSE259284)
RawGenotypes <- read.table(paste(path, "/GSM8113080_Crispr_Mouse_Batch1.mtDNA.RawGenotypes.Total.StrandBalance.txt.gz", sep = ""))
cov <- readRDS(paste0(path, '/GSM8113080_Crispr_Mouse_Batch1.mtDNA.depth.rds'))

# Fix CBs in coverage list
Dic_ATAC <- make_dic_ATAC()
cov$Total[[2]]$V1 <- to_ATAC(cov$Total[[2]]$V1, Dic_ATAC = Dic_ATAC)
cov <- list(Pos.MeanCov = cov$Total[[1]], Cell.MeanCov = cov$Total[[2]])

# Make VariantsGTSummary as in redeemR.read.trim()
GiveName <- c("UMI", "Cell", "Pos", "Variants", "Call", 
              "Ref", "FamSize", "GT_Cts", "CSS", "DB_Cts", "SG_Cts", 
              "Plus", "Minus", "Depth")
colnames(RawGenotypes) <- GiveName
RawGenotypes.annotated <- RawGenotypes %>% add_freq_raw() %>% make_position_df_3.5()
RawGenotypes.trimed <- filter(RawGenotypes.annotated, edge_dist >= edge_trim)
VariantsGTSummary.trimed <- GTSummary(RawGenotypes.trimed)
VariantsGTSummary <- GTSummary(RawGenotypes)
VariantsGTSummary <- VariantsGTSummary[, c("Var1", "Freq")] %>% 
  rename(Freq_before_trim = Freq) %>% 
  merge(VariantsGTSummary.trimed, ., by = "Var1") %>% 
  mutate(depth = depth - (Freq_before_trim - Freq))
VariantsGTSummary$hetero <- with(VariantsGTSummary, Freq/depth)

# Add attributes
attr(VariantsGTSummary, "thr") <- thr
attr(VariantsGTSummary, "path") <- path
attr(VariantsGTSummary, "edge_trim") <- edge_trim
attr(VariantsGTSummary, "depth") <- cov

# Create filtered redeemR object from raw (trimmed) basecall table
redeemR <- Create_redeemR_model(VariantsGTSummary, qualifiedCellCut=10, VAFcut=1, Cellcut=2)

# Filter redeem calls, last redeemR utils to do that
redeemR <- clean_redeem(redeemR, fdr=0.05)
redeemR <- clean_redeem_removehot(redeemR)

# Re-fix barcodes names
basecalls <- redeemR@GTsummary.filtered
meta <- redeemR@CellMeta
variants <- redeemR@V.fitered

# To RNA
Dic_RNA <- make_dic_RNA()
meta$Cell <- paste0(to_RNA(meta$Cell, Dic_RNA = Dic_RNA, PostFix = T), '-1')
basecalls$Cell <- paste0(to_RNA(basecalls$Cell, Dic_RNA = Dic_RNA, PostFix = T), '-1')
basecalls$Var1 <- paste0(basecalls$Cell, basecalls$Variants)

# Write out
write.csv(cov$Cell.MeanCov, paste0(path, '/coverage.csv'))
write.csv(basecalls, paste0(path, '/filtered_basecalls.csv'))
write.csv(meta, paste0(path, '/cell_meta.csv'))
write.csv(variants,  paste0(path, '/var_meta.csv'))


##
