#! /usr/bin/env/ Rscript

##when import to linux. vi file. and then :set ff = unix, :wq to save

library("tidyverse")
library("data.table")
library("dplyr")
library("tidyr")


##this script is dedicated to analyse the methylation of repetitive elements. It uses the repeat database rmsk as downloaded from ucsc.
##In the commandline you can specify which repeatclass you want to analyse. This refers to the repeat family such as ERVK, ERV1 etc.


##function to calculate avg methylation of a loci
loci_only_methyl_average <- function(df, chrom, loci_start, loci_end) {
  #Get from the merged dataframe, the whole region to work with.
  loci_only <- df[start >= loci_start & end <= loci_end & chrom == chromosome, 4:length(df)]
  avg_methyl_loci <- colMeans(loci_only, na.rm = T)
  avg_methyl_loci <- unname(avg_methyl_loci)
  return(avg_methyl_loci)
}




calculate_methylation_loci <- function(annotation_df, df) {
  annotation_df <- annotation_df[,1:3]
  ##make vectors from the cpg_df for usage in a loop
  chromosome_vec <- annotation_df[,1]
  start_vec <- annotation_df[,2]
  end_vec <- annotation_df[,3]
  
  loci_methylation_df <- matrix(nrow = length(chromosome_vec) ,ncol = length(df) - 3, data = NA)
  for (i in 1:length(chromosome_vec)) {
    annotation_methylation <- loci_only_methyl_average(df = df, chrom = chromosome_vec[i],
                                                       loci_start = start_vec[i], loci_end = end_vec[i])
    loci_methylation_df[i,] <- annotation_methylation
  }
  
  colnames(loci_methylation_df) <- colnames(df[,4:length(df)])
  loci_methylation_df <- as.data.frame(loci_methylation_df)
  
  return(loci_methylation_df)
}

calculate_repeat_methylation <- function(df, repeat_name) {
  repeats <- readRDS("rmsk.rds")
  df <- readRDS(df)
  df <- as.data.table(df)
  specific_repeat <- repeats[repeats$repFamily == repeat_name,]
  unique_names_ervk <- unique(specific_repeat$repName)
  print(paste("Total repeats are ",length(unique_names_ervk), sep = ""))
  
  repeats_methylation_df <- matrix(nrow = length(unique_names_ervk) ,ncol = length(df) -2, data = NA)
  for (i in 1:length(unique_names_ervk)) {
    print(i)
    name_spec_ervk_df <- specific_repeat[specific_repeat$repName == unique_names_ervk[i],]
    rep_methyl <- calculate_methylation_loci(name_spec_ervk_df, df)
    mean_rep_methyl <- colMeans(rep_methyl, na.rm = T)
    mean_rep_methyl <- unname(mean_rep_methyl)
    mean_rep_methyl <- c(unique_names_ervk[i], mean_rep_methyl)
    repeats_methylation_df[i,] <- mean_rep_methyl

    
  }
  repeats_methylation_df <- as.data.frame(repeats_methylation_df)
  colnames(repeats_methylation_df) <- c("repeat", colnames(df[,4:length(df)]))
  return(repeats_methylation_df)
  
}
##the command line arguments specified
##usage of this script: Rscript repeat_methylation.R methylation_dataframe.rds repeat_class output.rds
##as a measure of how long it takes the program will output a count for the loci. Depending on how many loci you supply you can see how long it will take. Easy to save this output to a file by command > log.txt &
args = commandArgs(trailingOnly = TRUE)
dataframe_touse <- args[1]
repeatname_touse <- args[2]
output <- args[3]
repeats_methylation <- calculate_repeat_methylation(df = dataframe_touse, repeat_name = repeatname_touse)



saveRDS(repeats_methylation,
        file = output)

