#!/usr/bin/env Rscript

# by Remi Allio, Damien de Vienne and Theo Tricou

args = commandArgs(trailingOnly=TRUE)

# load required packages
library(ape)
library(phylter)

#load arguments
trees<-read.tree(args[1])
# print(trees)

# bvalue: If X is a list of trees, nodes with a support below bvalue will be
# collapsed prior to the outlier detection.
bvalue<-as.numeric(args[2])

# distance: If X is a list of trees, type of distance used to compute the
# pairwise matrices for each tree. Can be “patristic” (sum of branch lengths
# separating tips, the default) or “nodal” (number of nodes separating tips).
distance<-as.character(args[3])

# k: Strength of outlier detection. The higher this value the less outliers
# detected.
k<-as.numeric(args[4])

# k2: Same as k for complete gene outlier detection. To preserve complete genes
# from being discarded, k2 can be increased. By default, k2 = k.
k2<-as.numeric(args[5])

# Norm: Should the matrices be normalized prior to the complete analysis and
# how. If “median”, matrices are divided by their median; if “mean”, they are
# divided by their mean; if “none”, no normalization if performed. Normalizing
# ensures that fast-evolving (and slow-evolving) genes are not treated as
# outliers. Normalization by median is a better choice as it is less sensitive
# to outlier values.
norm<-as.character(args[6])

# Norm.cutoff: Value of the median (if Norm = "median") or the mean (if Norm =
# "mean") below which matrices are simply discarded from the analysis. This
# prevents dividing by 0, and allows getting rid of genes that contain mostly
# branches of length 0 and are therefore uninformative anyway. Discarded genes,
# if any, are listed in the output (out$DiscardedGenes).
normcutoff<-as.numeric(args[7])

# gene.names: List of gene names used to rename elements in X. If NULL (the
# default), elements are named 1,2,…,length(X).
names<-args[8]

# test.island: If TRUE (the default), only the highest value in an island of
# outliers is considered an outlier. This prevents non-outliers hitchhiked by
# outliers to be considered outliers themselves.
island<-as.logical(args[9])

# verbose: If TRUE (the default), messages are written during the filtering
# process to get information on what is happening.
verbose<-as.logical(args[10])

# stop.criteria: The optimization stops when the gain (quality of compromise)
# between round n and round n+1 is smaller than this value. Default to 1e-5.
stopc<-as.numeric(args[11])

# InitialOnly: Logical. If TRUE, only the Initial state of the data is computed.
InitialOnly<-as.logical(args[12])

# normalizeby: Should the gene x species matrix be normalized prior to outlier
# detection, and how.
normby<-as.character(args[13])

# parallel: Logical. Should the computations be parallelized when possible?
# Default to TRUE. Note that the number of threads cannot be set by the user
# when parallel = TRUE. It uses all available cores on the machine.
parallel_task<-as.logical(args[14])

# pdfreport.file: If report=TRUE, name of the pdf file where the report is
# written. Default to report.pdf
report_to_pdf<-as.logical(args[15])

# Job name for outputs.
output_job<-as.character(args[16])

cat(paste("phylter(trees, bvalue = ", bvalue, ", distance = ", distance,
  ", k = ", k, ",k2 = ", k2, ", Norm = ", norm, ", Norm.cutoff = ", normcutoff,
  ", gene.names = ", names, ", test.island = ", island, ", verbose = ", verbose,
  ", stop.criteria = ", stopc, ", InitialOnly = ", InitialOnly,
  ", normalizeby = ", normby, ", parallel = ", parallel_task, ")\n", sep=""))

results<-phylter(trees, bvalue = bvalue, distance = distance, k = k, k2 = k2,
  Norm = norm, Norm.cutoff = normcutoff, gene.names = names,
  test.island = island, verbose = verbose, stop.criteria = stopc,
  InitialOnly = InitialOnly, normalizeby = normby, parallel = parallel_task)

write.phylter(results, file=paste(output_job, "/phylter.out", sep=""),
  pdfreport=report_to_pdf, pdfreport.file=paste(output_job, "/report.pdf", sep=""))
