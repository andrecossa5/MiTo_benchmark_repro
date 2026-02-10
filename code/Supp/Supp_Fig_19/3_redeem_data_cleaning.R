# Prep redeem mouse data
library(tidyverse)
library(redeemR)


##


# Utils
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


filter_v4_pipeline <- function(RawGenotypes, coverage.table, edge_trim=5, thr="Total") {
  
  # Trim RawGenotypes and create VariantsGTSummary
  RawGenotypes.annotated <- RawGenotypes %>% add_freq_raw() %>% make_position_df_3.5()
  RawGenotypes.trimed <- filter(RawGenotypes.annotated, edge_dist >= edge_trim)
  VariantsGTSummary.trimed <- GTSummary(RawGenotypes.trimed)
  VariantsGTSummary <- GTSummary(RawGenotypes)
  VariantsGTSummary <- VariantsGTSummary[, c("Var1", "Freq")] %>% 
    rename(Freq_before_trim = Freq) %>% 
    merge(VariantsGTSummary.trimed, ., by = "Var1") %>% 
    mutate(depth = depth - (Freq_before_trim - Freq))
  VariantsGTSummary$hetero <- with(VariantsGTSummary, Freq/depth)
  attr(VariantsGTSummary, "thr") <- thr
  attr(VariantsGTSummary, "edge_trim") <- edge_trim
  attr(VariantsGTSummary, "depth") <- coverage.table
  
  # Create redeemR object
  redeemR <- Create_redeemR_model(VariantsGTSummary, qualifiedCellCut=10, VAFcut=1, Cellcut=2)
  
  # Filter calls with the last redeemR utils
  redeemR <- clean_redeem(redeemR, fdr=0.05)
  redeemR <- clean_redeem_removehot(redeemR)
  
  # Format counts table
  counts <- redeemR@GTsummary.filtered %>% select(Cell, Variants, Freq, depth) %>% rename(AD_trimmed=Freq, DP=depth)
  counts$Variants <- str_replace(counts$Variants, "_([ACGTN])_([ACGTN])$", "_\\1>\\2")
  
  return(counts)
  
}


##


# Read data

# Path
path <- '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro/data/redeem_mouse/'

# N.B: The code below substitute redeemR.read.trim(), which cannot
# handle the Cell and Site coverage data provided by the authors in GSE259284 (RDS files, with RNA cell names)

# Read files (GEO: GSE259284)
RawGenotypes <- read.table(paste(path, "/GSE259284_RAW/GSM8113080_Crispr_Mouse_Batch1.mtDNA.RawGenotypes.Total.StrandBalance.txt.gz", sep = ""))
cov <- readRDS(paste0(path, '/GSE259284_RAW/GSM8113080_Crispr_Mouse_Batch1.mtDNA.depth.rds'))
cas9_table <- read.table(paste(path, "/redeem_mouse_batch1_allele_table.tsv", sep = ""))

# Fix MT-library CBs from ATAC --> RNA
GiveName <- c("UMI", "Cell", "Pos", "Variants", "Call", 
              "Ref", "FamSize", "GT_Cts", "CSS", "DB_Cts", "SG_Cts", 
              "Plus", "Minus", "Depth")
colnames(RawGenotypes) <- GiveName
Dic_RNA <- make_dic_RNA()
RawGenotypes$Cell <- to_RNA(RawGenotypes$Cell, Dic_RNA)

# Fix MT-coverage cell names (already RNA: remove -1)
cov$Total[[2]]$V1 <- sub('-1', '', cov$Total[[2]]$V1)

# Demultiplex RawGenotypes.Total.StrandBalance and coverage tables
# Trim RawGenotypes.Total.StrandBalance fragment-end basecalls
# Filter RawGenotypes with redeemR V4 filter
# Store RawGenotypes, coverage tables, and V4 filtered AD, DP counts for MiTo AFM building 

for (x in unique(cas9_table$sample)) {
  
  # Subset RawGenotypes and meanCellCov tables
  cells <- cas9_table[cas9_table$sample == x, "cellBC"]
  genotypes <- RawGenotypes %>% filter(Cell %in% cells)
  meanSiteCov <- cov$Total[[1]]
  meanCellCov <- cov$Total[[2]] %>% filter(V1 %in% cells)

  # Filter MT-SNVs (redeemR v4 filter) and return filtered AD/DP counts
  filtered_counts <- filter_v4_pipeline(genotypes, 
                                        edge_trim=0,
                                        list(Pos.MeanCov=meanSiteCov, Cell.MeanCov=meanCellCov))
  filtered_counts.trimmed <- filter_v4_pipeline(genotypes, 
                                        edge_trim=5,
                                        list(Pos.MeanCov=meanSiteCov, Cell.MeanCov=meanCellCov))
  # Write out everything
  write.table(genotypes, 
              paste0(path, x, '/RawGenotypes.Total.StrandBalance'), 
              col.names=FALSE, row.names=FALSE, sep='\t') 
  write_tsv(filtered_counts, paste0(path, x, '/FilteredCounts')) 
  write_tsv(filtered_counts.trimmed, paste0(path, x, '/FilteredCounts.trimmed')) 
  write_tsv(meanSiteCov %>% rename(pos=V2, coverage=meanCov), paste0(path, x, '/meanSiteCov')) 
  write_tsv(meanCellCov %>% rename(cell=V1, coverage=meanCov), paste0(path, x, '/meanCellCov')) 
  
}


##
