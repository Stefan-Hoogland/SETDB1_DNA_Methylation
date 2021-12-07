#! /usr/bin/env/ Rscript

library("tidyverse")
library("data.table")
library("dplyr")
library("tidyr")


##The below two functions work together to take a bed file and the methylation dataframe. Using this bed file it will calculate the average methylation for each loci.
##It will also give the mean methylation of the surrounding areas -10kb to +10kb in increments of 500bp


##This function aids the below function. It runs on the background.
loci_methyl_average <- function(df, chrom, loci_start, loci_end) {
  #Get from the merged dataframe, the whole region to work with.
  methyl_df <- df[start >= loci_start - 10000 & end <= loci_end + 10000 & chromosome == chrom,]
  
  ##Get a vector of only the loci methylation
  loci_only <- methyl_df[start >= loci_start & end <= loci_end, 4:16]
  avg_methyl_loci <- colMeans(loci_only, na.rm = T)
  avg_methyl_loci <- unname(avg_methyl_loci)
  
  ##specify a range
  range <- seq(0,10000, by = 500)
  ##Create an empty dataframe where all the mean methylations of each range can be put
  left_mean_ranges_df <- data.frame(range_name = character(length(range)-1), serum_control = integer(length(range)-1), serum_day4 = integer(length(range)-1), 
                                    serum_day6 = integer(length(range)-1), serum_day8 = integer(length(range)-1), serum_day10 = integer(length(range)-1),
                                    t2i_control = integer(length(range)-1), t2i_day1 = integer(length(range)-1),
                                    t2i_day2 = integer(length(range)-1), t2i_day4 = integer(length(range)-1),
                                    serumto2i_day1 = integer(length(range) -1), serumto2i_day2 = integer(length(range)-1),
                                    serumto2i_day3 = integer(length(range)-1), serumto2i_day4 = integer(length(range)-1))
  
  right_mean_ranges_df <- data.frame(range_name = character(length(range)-1), serum_control = integer(length(range)-1), serum_day4 = integer(length(range)-1), 
                                     serum_day6 = integer(length(range)-1), serum_day8 = integer(length(range)-1), serum_day10 = integer(length(range)-1),
                                     t2i_control = integer(length(range)-1), t2i_day1 = integer(length(range)-1),
                                     t2i_day2 = integer(length(range)-1), t2i_day4 = integer(length(range)-1),
                                     serumto2i_day1 = integer(length(range) -1), serumto2i_day2 = integer(length(range)-1),
                                     serumto2i_day3 = integer(length(range)-1), serumto2i_day4 = integer(length(range)-1))
  for (i in 1:length(range)) {
    #failsave for when the last i gives an error
    if (i+1 > length(range)) {
      break
    }
    else {
      #left side of the loci
      left_ranges_df <- methyl_df[methyl_df$start >= loci_start - range[i + 1] & methyl_df$end < loci_start - range[i],4:16]
      left_mean_range_vec <- unname(colMeans(left_ranges_df, na.rm = T))
      left_range_name <- paste("-",range[i],":","-",range[i+1], sep = '')
      left_mean_ranges_df[i,1] = left_range_name
      left_mean_ranges_df[i,2:14] = left_mean_range_vec
      
      #loci itself
      cpg_specific_df <- methyl_df[methyl_df$start >= loci_start & methyl_df$end <= loci_end, 4:16]
      mean_cpg_df <- data.frame(colMeans(cpg_specific_df, na.rm = T))
      mean_cpg_df <- as.data.frame(t(mean_cpg_df))
      mean_cpg_df$range_name = "loci"
      
      #right side of the loci
      right_ranges_df <- methyl_df[methyl_df$start >= loci_end + range[i] & methyl_df$end < loci_end + range[i +1],4:16]
      right_mean_range_vec <- unname(colMeans(right_ranges_df, na.rm = T))
      right_range_name <- paste(range[i],":",range[i+1], sep = '')
      right_mean_ranges_df[i,1] = right_range_name
      right_mean_ranges_df[i,2:14] = right_mean_range_vec
    }
  }
  ranges_df <- bind_rows(left_mean_ranges_df, mean_cpg_df, right_mean_ranges_df)
  rownames(ranges_df) <- NULL
  methyl_list <- list(ranges_df, avg_methyl_loci)
  
  
  
  return(methyl_list)
}






##This function will perform the calculation. The end result is a dataframe with for each loci the mean methylation and surrounding methylation in increments of 500bp.
##file_location refers to the location of the bed_file. It is important that the first three columns are: Chromosome, start, end.
##df refers to the methylation dataframe created using bisulfite sequencing data. it is stored as a RDS file for quick read in by R
calculate_methylation <- function(file_location, df) {
  
  ##load in the dataframe from the file_location
  annotation_df <- read.table(file_location, sep = "\t", header = F, stringsAsFactors = F)
  annotation_df <- na.omit(annotation_df)
  print(head(annotation_df))
  df <- readRDS(df)

  
  ##make vectors from the cpg_df for usage in a loop
  chromosome_vec <- annotation_df[,1]
  start_vec <- as.integer(annotation_df[,2])
  end_vec <- as.integer(annotation_df[,3])
  
  ##convert dataframe to data.table
  dt <- data.table(df)
  
  ##create a list for the dataframes of each annotation using the function above
  annotation_ranges_list <- list()
  loci_methylation_df <- matrix(nrow = length(chromosome_vec) ,ncol = 13, data = NA)
  for (i in 1:length(chromosome_vec)) {
    print(i)
    annotation_methyl_list <- loci_methyl_average(df = dt, chrom = chromosome_vec[i],
                                                  loci_start = start_vec[i], loci_end = end_vec[i])
    annotation_methyl_ranges <- annotation_methyl_list[[1]]
    annotation_ranges_list[[i]] <- annotation_methyl_ranges
    loci_methylation <- annotation_methyl_list[[2]]
    loci_methylation_df[i,] <- loci_methylation
  }
  
  colnames(loci_methylation_df) <- colnames(dt[,4:16])
  loci_methylation_df <- as.data.frame(loci_methylation_df)
  
  ##calculate the mean methylation of each cell from each dataframe in the list
  all_annotation_methyl_range_df <- rbindlist(annotation_ranges_list)[,lapply(.SD, mean, na.rm = T), list(range_name)]
  
  
  ##rename the range_name for usage in plot, and order the dataframe
  range_name <- c(-500,-1000,-1500,-2000,-2500,-3000,-3500,-4000,-4500,-5000,-5500,-6000,-6500,-7000,-7500,-8000,-8500,-9000,
                  -9500,-10000,0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000)
  all_annotation_methyl_range_df$range_name <- range_name
  all_annotation_methyl_range_df <- all_annotation_methyl_range_df[order(range_name),]
  final_methylation_list <- list(all_annotation_methyl_range_df, loci_methylation_df)
  
  return(final_methylation_list)
}

##the command line arguments specified
##usage of this script: Rscript annotation_methylation.R file.bed methylation_dataframe.rds output.rds
##as a measure of how long it takes the program will output a count for the loci. Depending on how many loci you supply you can see how long it will take. Easy to save this output to a file by command > log.txt &
args = commandArgs(trailingOnly = TRUE)
annotation_file <- args[1]
dataframe_touse <- args[2]
output <- args[3]
annotation_methylation <- calculate_methylation(annotation_file, df = dataframe_touse)



saveRDS(annotation_methylation,
        file = output)
