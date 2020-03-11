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
# new_names = unlist(lapply(tr$node.label, function(x)strsplit(x, "_")[[1]][1]))
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

# data = read.table("AB_BA.txt", sep = "\t", h = T)
# internal = read.table("internal_ABBA_BABA.txt", sep = "\t", h = T)







###### version of the code that remove snp shared by an ancestor in all its leaves.


tr = read.tree('spe_rename')
plot(tr, show.node=T)
data = read.table('data',h = T,row.names = 1, sep = "\t")

aa=c()
for (i in 1:ncol(data)) {
  if (length(unique(data[,i])) == 1 ) {
    aa = c(aa,i)
  }
}
data = data[,-aa]
dim(data)


liste_node = c(tr$tip.label, tr$node.label)


for (i in length(tr$node.label):2) {
  name = tr$node.label[i]
  data[name,] <- rep(NA, ncol(data))
  leaf_id = Descendants(tr, which(liste_node == name), type = "children")
  leaf = liste_node[leaf_id]
  sub = data[leaf,]
  for (i in 1:ncol(sub)) {
    if (length(unique(sub[,i])) == 1) {
      data[name, i] = unique(sub[,1])
      data[leaf, i] = NA
    }
  }
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
  print(length(arr))
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




topology = getquatuors_internal(tr)
nrow(topology)


is_aaaa <- function(seg_site) {if(!NA %in% seg_site){if(sum(as.numeric(seg_site)) == 4 || sum(as.numeric(seg_site)) == 0){return(1)} else {return(0)}} else {return(0)}}
is_baaa <- function(seg_site) {if(!NA %in% seg_site){if(sum(as.numeric(seg_site)) == 1) {return(1)} else {return(0)}}else {return(0)}}
is_abba <- function(seg_site) {if(!NA %in% seg_site){if(sum(as.numeric(seg_site)) == 2 & seg_site[1] == seg_site[4]){return(1)} else {return(0)}} else {return(0)}}
is_baba <- function(seg_site) {if(!NA %in% seg_site){if(sum(as.numeric(seg_site)) == 2 & seg_site[1] == seg_site[3]){return(1)} else {return(0)}}else {return(0)}}
is_bbaa <- function(seg_site) {if(!NA %in% seg_site){if(sum(as.numeric(seg_site)) == 2 & seg_site[1] == seg_site[2]){return(1)} else {return(0)}}else {return(0)}}


D_stat <- function(data, quatuor){
  temp = data[quatuor,]
  aaaa <- sum(apply(temp, 2, function(x) is_aaaa(x)))
  baaa <- sum(apply(temp, 2, function(x) is_baaa(x)))
  bbaa <- sum(apply(temp, 2, function(x) is_bbaa(x)))
  abba <- sum(apply(temp, 2, function(x) is_abba(x)))
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


quatuor = c("maga1260", "mait1250", "bhoj1244", "norw1258")


quatuor = c("maga1260","mait1250","bhoj1244","N4_jaun1243_kuma1273_nepa1254")
temp = data[quatuor,]

quatuor= topology[1,]
seg_site = temp[,52]
is_aaaa(seg_site)

is_aaaa <- function(seg_site) {if(!NA %in% seg_site){if(sum(as.numeric(seg_site)) == 4 || sum(as.numeric(seg_site)) == 0){return(1)} else {return(0)}} else {return(0)}}


quatuor = c("maga1260","mait1250","bhoj1244","N4_jaun1243_kuma1273_nepa1254")
temp = data[quatuor,]

for (i in 1:ncol(temp)) {
  seg_site = temp[,i]
  print(i)
  print(is_abba(seg_site))
  print(is_baaa(seg_site))
}
seg_site = temp[,52]
print(is_abba(seg_site))
print(is_baaa(seg_site))



D_stat(data, quatuor)


results <- as.data.frame(t(apply(topology, 1, function(x) D_stat(data, x))))
write.table(results, "ancestor_internal_ABBA_BABA.txt", sep = "\t", row.names = F, append = F, quote = F)


signi = results[as.numeric(as.character( results$Pvalue)) < 0.05,]
write.table(signi, "signi_ancestor_internal_ABBA_BABA.txt", sep = "\t", row.names = F, append = F, quote = F)


signi[order(signi$P3),]


N11

quatuor = c( "hind1269", "awad1243", "N11_panj1256_sind1272","norw1258")
temp = data[quatuor,]


for (i in 1:ncol(temp)) {
  seg_site = temp[,i]
  print(i)
  print(is_abba(seg_site))
  print(is_baba(seg_site))
}
