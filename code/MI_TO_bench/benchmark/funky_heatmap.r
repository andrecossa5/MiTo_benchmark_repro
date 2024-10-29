# Funky heatmap

library(tidyverse)
library(funkyheatmap)

# data("dynbenchmark_data")

my_table <- read.csv('/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/results/MI_TO_bench/phylo_inference/grouped.csv', row.names = 1)
my_table <- my_table %>% add_column(id = c(1:5,66:70) %>% as.character, .before = 1)

col.names <- c(
  "", "Preprocessing", "Binarization", "Min confident AF", "Min n ALT", "Min n +cells", 
  "Mut Quality", "GBC", "Tree structure", "Connectedness", "Variation", "Yield", "Overall",
  "n dbSNP", "n REDIdb", "Mut signature", 
  "Clonal-biased MT-SNVs", "AUPRC", "ARI", "NMI", 
  "Tree vs MT-SNVs corr", 
  "Density", "Transitivity", "Mean path length", "Average degree", "LLC", 
  "Aplotype redundancy", 'n MT-SNVs',
  "Median n MT-SNVs (cell)", "n GBC clones", "n cells"
)

geoms <- c(rep("text",6), rep("bar", 7), rep("funkyrect", 18))

palettes <- c(
  rep(NA,6), 
  "Mutation Quality", "Association with GBC", "Tree Structure","Connectedness", "Var", "Yield", "Overall", 
  rep( "Mutation Quality", 3),
  rep("Association with GBC",4),
  "Tree Structure",
  rep("Connectedness", 5), 
  rep("Var",2),
  rep("Yield", 3)
)

col.groups <- c(
  rep("Parameters",6), rep("Summary",7), rep("Mutation Quality", 3), 
  rep("Association with GBC", 4), "Tree Structure", rep("Connectedness", 5), rep("Var",2), rep("Yield", 3)
)

column_info <- tibble(
  group=col.groups,
  id=colnames(my_table), 
  name=col.names,
  geom=geoms,
  palette=palettes,
  legend=rep(FALSE, 31),
  width=c(c(2,5,3.5,1,1,1), rep(1.5,7), rep(1,18)),
  overlay=rep(FALSE,31)
)

yellow_gradient <- c(
  "#333300", "#4D4D00", "#666600", "#808000", "#999900", 
  "#B3B300", "#CCCC00", "#E6E600", "#FFFF00", "#FFFF33", 
  "#FFFF66", "#FFFF99", "#FFFFCC", "#FFFFE6", "#FFFFF0"
)

purple_gradient <- c(
  "#2D003D", "#400055", "#52006E", "#660088", "#7900A2", 
  "#8C00BB", "#A100D4", "#B400ED", "#C766FF", "#D599FF", 
  "#E3CCFF", "#F2E6FF", "#F8F0FF", "#FAF5FF", "#FCFBFF"
)

brown_gradient <- c(
  "#331900", "#4D2600", "#663300", "#804000", "#994D00", 
  "#B35A00", "#CC6600", "#E67300", "#FF8000", "#FF9933", 
  "#FFB266", "#FFCC99", "#FFE6CC", "#FFF0E6", "#FFF7F0"
)

palettes <- list(
  "Mutation Quality" = "Blues",
  "Association with GBC" = "Reds",
  "Tree Structure" = "Greens",
  "Connectedness" = brown_gradient, 
  "Var" = yellow_gradient, 
  "Yield" = purple_gradient,
  "Overall" = "Greys"
)

col.groups.df <- tibble(
  group=col.groups %>% unique,
  palette=c("Overall", "Overall", "Mutation Quality", "Association with GBC", "Tree Structure", "Connectedness", "Var", "Yield"),
  level1=col.groups %>% unique
)


##


g <- funkyheatmap::funky_heatmap(
  data = my_table,
  column_info = column_info,
  column_groups = col.groups.df,
  palettes = palettes,
  position_args = position_arguments(col_bigspace = 1, col_annot_offset = 3, col_annot_angle = 45, expand_ymax = 10, expand_xmax = 5, expand_xmin=5, expand_ymin=10),
  scale_column = FALSE  # This will scale numerical columns between their min and max
)

ggsave('/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/results/MI_TO_bench/benchmark/heatmap_MDA_lung.pdf', g, device = cairo_pdf, width = g$width, height = g$height)


##
