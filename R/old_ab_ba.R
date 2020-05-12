t = tr('spe_tree', format = 1 )


count=1
for i in t.traverse():
  if not i.is_leaf():
    ll=["N", str(count)]
    count += 1
    for j in  i.get_leaves() :
      ll.append(j.name)
    i.name = "_".join(ll)



for i in t.traverse(): i.name


R

require('ape')
tr = read.tree('spe_tree')

validateandreorder<-function(arr, dist) {
  # created by Damien
  submat<-dist[arr,arr]
  if (sum(submat==max(submat))==6 & length(unique(levels( as.factor(submat)))) == 4) {
    diag(submat)<-1
    return(names(sort(apply(submat,1,prod))))
  }
}

getquatuors<-function(tr) {
  # created by Damien, modified by Theo
  # dist<-cophenetic(compute.brlen(tr)) # compute.brlen used to get rid of uneven branch length
  dist<-round(cophenetic(tr),3) # this is a dangerous thing to do
  allquat<-combn(tr$tip.label,4)
  RES<-do.call(rbind,apply(allquat, 2, function(x) validateandreorder(x, dist)))
  return(RES)
}



data = read.table('data',sep = '\t', h = T)




allquat<-combn(tr$tip.label,4)
dist<-cophenetic(tr) # this is a dangerous thing to do
dist<-round(cophenetic(tr),3) # this is a dangerous thing to do

arr = allquat[,1]
(submat<-dist[arr,arr])


arr = allquat[,20475]

(submat<-dist[arr,arr])
if (sum(submat==max(submat))==6) {

topology = getquatuors(tr)



data = read.table('data',sep = '\t', h = T)

is_aaaa <- function(seg_site) {if(sum(seg_site) == 4 || (sum(seg_site) == 0)) {return(1)} else {return(0)}}
is_baaa <- function(seg_site) {if(sum(seg_site) == 1) {return(1)} else {return(0)}}
is_abba <- function(seg_site) {if(sum(seg_site) == 2 & seg_site[1] == seg_site[4]){return(1)} else {return(0)}}
is_baba <- function(seg_site) {if(sum(seg_site) == 2 & seg_site[1] == seg_site[3]){return(1)} else {return(0)}}
is_bbaa <- function(seg_site) {if(sum(seg_site) == 2 & seg_site[1] == seg_site[2]){return(1)} else {return(0)}}


D_stat <- function(stat_simulation, quatuor){
  q = topology[1,]
  p1 = quatuor[1]
  p2 = quatuor[2]
  p3 = quatuor[3]
  p4 = quatuor[4]
  temp = as.data.frame(rbind(data[data$Glottocode == c(p1),],
  data[data$Glottocode == c(p2),],
  data[data$Glottocode == c(p3),],
  data[data$Glottocode == c(p4),]))
  aaaa <- sum(apply(temp[,-1], 2, function(x) is_aaaa(x)))
  baaa <- sum(apply(temp[,-1], 2, function(x) is_baaa(x)))
  bbaa <- sum(apply(temp[,-1], 2, function(x) is_abba(x)))
  abba <- sum(apply(temp[,-1], 2, function(x) is_bbaa(x)))
  baba <- sum(apply(temp[,-1], 2, function(x) is_baba(x)))
  if ((abba + baba) != 0) {D = (abba - baba) / (abba + baba)} else {D = 0}
  data <- c("P1" = quatuor[1], "P2" = quatuor[2], "P3" = quatuor[3], "P4" = quatuor[4],
  "aaaa" = aaaa, "baaa" = baaa, "bbaa" = bbaa, "abba" = abba, "baba" = baba, "Dstat" = D,
  "Pvalue" = binom.test(c(abba, baba), p = 0.5, conf.level= 0.95)$p.value)
  return(data)
}




results <- as.data.frame(t(apply(topology, 1, function(x) D_stat(data, x))))
write.table(results, "AB-BA.txt", sep = "\t", row.names = F, append = F, quote = F)



data = read.table("AB-BA.txt", sep = "\t", h = T)

df = data[data$Pvalue < 0.01,]

df[order(df$P3, df$P2),]
