#!/usr/bin/env Rscript
# Théo

require('ape')
require('phangorn')

old_tongue<-function(name, data) {
  names = c(strsplit(name, "_")[[1]][-1])
  sub = data[names,]
  res <- apply(sub, 2, function(x) sound(x))
  return(res)
}

sound <- function(col) {
  if (length(unique(col)) == 1) {
    return(unique(col))
  } else {return(NA)}
}

getquatuors_internal <- function(tr) {
  dist<-round(dist.nodes(tr),2)
  order = c(tr$tip.label, tr$node.label)
  allquat<-combn(order[-which(order == "World")],4)
  rownames(dist) = order
  colnames(dist) = order
  RES<-do.call(rbind,apply(allquat, 2, function(x) validateandreorder(x, dist, tr)))
  return(RES)
}

validateandreorder<-function(arr, dist, tr) {
  if (length(which(startsWith(arr, "N")) == TRUE) != 0) {
    test = which(order %in% arr) %in% unlist(Descendants(tr, which(order %in% arr[startsWith(arr, "N")]), type = "all"))
  } else {
    test = FALSE
  }
  if (!TRUE %in% test) {
    array_leaves = sapply(arr, function(x) one_leaf(tr, x))
    submat<-dist[array_leaves,array_leaves]
    diag(submat)<-1
    if (length(unique(apply(submat,1,prod))) == 3) {
      or = names(sort(apply(submat,1,prod)))
      res = c(names(which(array_leaves == or[1])),names(which(array_leaves == or[2])),names(which(array_leaves == or[3])),names(which(array_leaves == or[4])))
      return(res)
    }
  }
}

one_leaf <-function(tree, name) {
  if (startsWith(name, "N")) {
    return (extract.clade(tree, which(order == name))$tip.label[1])
  } else {
    return(name)
  }
}

tr = read.tree('spe_rename')
data = read.table('data',h = T,row.names = 1, sep = "\t")

new_data =  t(rbind(sapply(tr$node.label, function(x) old_tongue(x, data))))
new_names = unlist(lapply(tr$node.label, function(x)strsplit(x, "_")[[1]][1]))
internal = rbind(data, new_data)

topology = getquatuors_internal(tr)
nrow(topology)

is_aaaa <- function(seg_site) {
  if(!NA %in% seg_site){
    if(sum(seg_site) == 4 || sum(seg_site) == 0){
      return(1)
    } else {
      return(0)
    }
  } else {
    return(0)
  }
}

is_baaa <- function(seg_site) {if(!NA %in% seg_site){if(sum(seg_site) == 1) {return(1)} else {return(0)}}else {return(0)}}
is_abba <- function(seg_site) {if(!NA %in% seg_site){if(sum(seg_site) == 2 & seg_site[1] == seg_site[4]){return(1)} else {return(0)}} else {return(0)}}
is_baba <- function(seg_site) {if(!NA %in% seg_site){if(sum(seg_site) == 2 & seg_site[1] == seg_site[3]){return(1)} else {return(0)}}else {return(0)}}
is_bbaa <- function(seg_site) {if(!NA %in% seg_site){if(sum(seg_site) == 2 & seg_site[1] == seg_site[2]){return(1)} else {return(0)}}else {return(0)}}


D_stat <- function(data, quatuor){
  temp = data[quatuor,]
  aaaa <- sum(apply(temp, 2, function(x) is_aaaa(x)))
  baaa <- sum(apply(temp, 2, function(x) is_baaa(x)))
  bbaa <- sum(apply(temp, 2, function(x) is_abba(x)))
  abba <- sum(apply(temp, 2, function(x) is_bbaa(x)))
  baba <- sum(apply(temp, 2, function(x) is_baba(x)))
  if ((abba + baba) != 0) {D = (abba - baba) / (abba + baba)} else {D = 0}
  data <- c("P1" = strsplit(quatuor[1], "_")[[1]][1],
  "P2" = strsplit(quatuor[2], "_")[[1]][1],
  "P3" = strsplit(quatuor[3], "_")[[1]][1],
  "P4" = strsplit(quatuor[4], "_")[[1]][1],
  "aaaa" = aaaa, "baaa" = baaa, "bbaa" = bbaa, "abba" = abba, "baba" = baba, "Dstat" = D,
  "Pvalue" = binom.test(c(abba, baba), p = 0.5, conf.level= 0.95)$p.value)
  return(data)
}


results <- as.data.frame(t(apply(topology, 1, function(x) D_stat(internal, x))))
write.table(results, "internal_ABBA_BABA.txt", sep = "\t", row.names = F, append = F, quote = F)

data = read.table("AB_BA.txt", sep = "\t", h = T)
internal = read.table("internal_ABBA_BABA.txt", sep = "\t", h = T)
