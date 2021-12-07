#! /usr/bin/env/ Rscript

library("tidyverse")
library("data.table")
library("dplyr")
library("tidyr")

##This script is very similar to annoation_methylation.R. However, it now only returns the mean methylation of the loci itself. This makes it much faster if you have no interest in the surrounding methylation.

##Function that takes dataframe, chr en loci start and end and gives average methylation within a specified range
loci_only_average <- function(df, chrom, loci_start, loci_end) {
  
  ##Get a vector of only the loci methylation
  loci_only <- df[start >= loci_start & end <= loci_end & chromosome == chrom, 4:16]
  avg_methyl_loci <- colMeans(loci_only, na.rm = T)
  avg_methyl_loci <- unname(avg_methyl_loci)
  
  return(avg_methyl_loci)
}


##function that accepts a annotation dataframe with the genomic locations in the first three columns
calculate_methylation <- function(annotation_df, df) {
  annotation_df <- read.table(annotation_df, header = F, sep = "\t", stringsAsFactors = F)
  annotation_df <- na.omit(annotation_df)
  print(head(annotation_df))
  
  ##make vectors from the cpg_df for usage in a loop
  chromosome_vec <- annotation_df[,1]
  start_vec <- as.integer(annotation_df[,2])
  end_vec <- as.integer(annotation_df[,3])
  
  ##convert dataframe to data.table
  df <- readRDS(df)
  dt <- data.table(df)
  
  ##calculate methylation of each loci
  loci_methylation_df <- matrix(nrow = length(chromosome_vec) ,ncol = 13, data = NA)
  for (i in 1:length(chromosome_vec)) {
    print(i)
    annotation_methyl <- loci_only_average(df = dt, chrom = chromosome_vec[i],
                                           loci_start = start_vec[i], loci_end = end_vec[i])
    loci_methylation_df[i,] <- annotation_methyl
  }
  ##make final output dataframe
  colnames(loci_methylation_df) <- colnames(df[,4:16])
  loci_methylation_df <- as.data.frame(loci_methylation_df)
  loci_methylation_df <- cbind(annotation_df, loci_methylation_df)
  
  
  return(loci_methylation_df)
}


##the command line arguments specified
##usage of this script: Rscript loci_only_methylation.R file.bed methylation_dataframe.rds output.rds
##as a measure of how long it takes the program will output a count for the loci. Depending on how many loci you supply you can see how long it will take. Easy to save this output to a file by command > log.txt &
args = commandArgs(trailingOnly = TRUE)
annotation_file <- args[1]
dataframe_touse <- args[2]
output <- args[3]
annotation_methylation <- calculate_methylation(annotation_file, df = dataframe_touse)



saveRDS(annotation_methylation,
        file = output)