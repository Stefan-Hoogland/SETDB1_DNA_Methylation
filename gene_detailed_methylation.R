#! /usr/bin/env/ Rscript

library("tidyverse")
library("data.table")
library("dplyr")
library("tidyr")

##This script takes the gene database from ucsc (refseq) and will calculate the mean methylation on regulatory region and gene body for each gene. note: to use this script you need to alter the file locations for the
##methylation dataframe and the genes dataframe.
##This script will take + and - strand into account.

genes <- read.table("genes_detailed.txt",
                    header = T, sep = "\t")##alter gene dataframe location as needed

print(head(genes))


##regulatory region = (first_exon_start -1000) till first_exon_end
##gene body = first exon end till tx end
gapped_merged <- readRDS("/scratch/shoogland/dna_methylation_data/new_analysis/merged_methylation/gapped_merged.rds")##alter methylation_dataframe location as needed
print(head(gapped_merged))
dt <- data.table(gapped_merged)


regulatory_region_methylation_df <- matrix(nrow = nrow(genes), ncol = 13, data = NA)
gene_body_methylation_df <- matrix(nrow = nrow(genes), ncol = 13, data = NA)
for (i in 1:nrow(genes)) {
  print(i)
  if (genes[i,5] == "+") {
    ##if gene is on positive strand
    regulatory_region <- dt[start >= genes[i,6] - 1000 & end <= genes[i,7] & chromosome == genes[i,2],4:16]
    gene_body <- dt[start >= genes[i,7] & end <= genes[i,4] & chromosome == genes[i,2],4:16]
  }
  if (genes[i,5] ==  "-") {
    ##if gene is on negative strand
    regulatory_region <- dt[start >= genes[i,8] & end <= genes[i,9] + 1000 & chromosome == genes[i,2],4:16]
    gene_body <- dt[start >= genes[i,3] & end <= genes[i,8] & chromosome == genes[i,2],4:16]
  }
  regulatory_region_methylation <- unname(colMeans(regulatory_region, na.rm = T))
  regulatory_region_methylation_df[i,] <- regulatory_region_methylation
  
  gene_body_methylation <- unname(colMeans(gene_body, na.rm = T))
  gene_body_methylation_df[i,] <- gene_body_methylation
}

regulatory_region_methylation_df <- data.frame(regulatory_region_methylation_df)
colnames(regulatory_region_methylation_df) <- c("serum_control", "serum_day4", "serum_day6", "serum_day8", "serum_day10", "t2i_control",
                                                "t2i_day1", "t2i_day2", "t2i_day4", "serumto2i_day1", "serumto2i_day2", "serumto2i_day3",
                                                "serumto2i_day4")
regulatory_region_methylation_df$name <- genes[,1]

gene_body_methylation_df <- data.frame(gene_body_methylation_df)
colnames(gene_body_methylation_df) <- c("serum_control", "serum_day4", "serum_day6", "serum_day8", "serum_day10", "t2i_control",
                                        "t2i_day1", "t2i_day2", "t2i_day4", "serumto2i_day1", "serumto2i_day2", "serumto2i_day3",
                                        "serumto2i_day4")
gene_body_methylation_df$name <- genes[,1]

saveRDS(regulatory_region_methylation_df, file = "regulatory_region_methylation.rds")
saveRDS(gene_body_methylation_df, file = "gene_body_methylation.rds")



