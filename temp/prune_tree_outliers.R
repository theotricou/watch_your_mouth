#!/usr/bin/env Rscript

# by  Damien de Vienne and Theo Tricou


args = commandArgs(trailingOnly=TRUE)

# load required packages
library(ape)

#load arguments

# Gene trees directory.
trees_dir<-as.character(args[1])

# Phylter output file containing outliers.
job_dir<-as.character(args[2])

#outlier data frame
outliers <- read.table(args[2], h=F)



if (nrow(outliers)>=1){
  # list gene tree files
  trees <- Sys.glob(paste(trees_dir,"/*", sep=''))
  # identify gene tree to be pruned
  to_prune <- unique(grep(paste(outliers[,1],collapse="|"), trees, value=TRUE))
  # identify gene tree taht stay the same
  same <- trees[!(trees %in% to_prune)]




  dir.create(file.path(paste(job_dir, "/trees_PhylteR", sep="")), showWarnings = FALSE)

  # this funciton
  prunetree <- function(tree_name){
    # this function use keep.tip from the ape package to prune tree with outliers
    gene_tree_name <- to_prune[grep(tree_name[1,1], to_prune)]
    tree <- read.tree(gene_tree_name)

    if (nrow(tree_name) != length(tree$tip.label)){
      to_remove <- tree_name[,2]
      pruned_tree <- keep.tip(tree, tree$tip.label[-which(tree$tip.label %in% to_remove)])
      if (length(pruned_tree$tip.label) > 2){
        write.tree(pruned_tree, file = paste(job_dir, "/trees_PhylteR/", tree_name[1,1], "_phylter.nwk", sep=""))
      }else{
        cat(paste("Tree from gene", tree_name[1,1], "was removed due to insufficient number of leaves remaining: nTip =", length(pruned_tree$tip.label, "\n")))
      }
    }else{
      cat(paste("Tree from gene", tree_name[1,1], "was removed due to insufficient number of leaves remaining: nTip = 0\n"))
    }
  }
  by(outliers, outliers[,1], function(x) prunetree(x))


  cat(paste("The remaining ", length(same) , " gene tree are unfiltered\n",  sep = ""))
  invisible(lapply(same, function(x) {
    cmd <- paste("cp ", x, " ", job_dir, "/trees_PhylteR/", sep ="")
    system(cmd)
    }
  ))
}else{
  cat("There is no outliers and no tree to prune\n")
}
