#!/usr/bin/env Rscript
# Theo

library("ape")

setwd("~/Downloads/6_work_in_progress/Glotto/4_taxons_test/")

file<-"ABIR"
file<-"GFEA"
file<-"PBIR"

flexi<-paste(file,"_lexicon.tsv",sep="")
fmorph<-paste(file,"_morphosyntax.csv",sep="")
fphono<-paste(file,"_phono.csv",sep="")
tree<-paste("tree_",file,".nex",sep="")

lexi<-t(read.table(paste(file,flexi,sep="/"), h=T,row.names = 1))
morph<-t(read.table(paste(file,fmorph,sep="/"),sep=",",h=T,row.names = 1))
phono<-t(read.table(paste(file,fphono,sep="/"),sep=",",h=T,row.names = 1))

t<-read.tree(paste(file,tree,sep="/"))
dist<-cophenetic(t)
quat<-combn(t$tip.label,4)
diag(dist)<-1
ord<-names(sort(apply(dist,1,prod)))

paternR<-function(data){
  if (length(unique(data)) == 1){return("aaaa")} # AAAA
  else{
    if (sum(data) == 2){
      if (data[ord[1]] == data[ord[2]]){return("bbaa")} # BBAA
      else if (data[ord[1]] == data[ord[3]]){return("baba")} # BABA
      else if (data[ord[1]] == data[ord[4]]){return("abba")} # ABBA
    }
    else if (length(unique(data[c(ord[1],ord[2],ord[3])])) == 1){return("bbba")} # BBBA
    else if (length(unique(data[c(ord[1],ord[2],ord[4])])) == 1){return("aaba")} # AABA
    else if (length(unique(data[c(ord[1],ord[3],ord[4])])) == 1){return("abaa")} # ABAA
    else if (length(unique(data[c(ord[2],ord[3],ord[4])])) == 1){return("baaa")} # BAAA
  }
}


lex<-apply(lexi,1,function(x) paternR(x))
mor<-apply(morph,1,function(x) paternR(x))
pho<-apply(phono,1,function(x) paternR(x))

tlex<-as.data.frame(lex[which(lex == "abba" | lex == "baba" )])
colnames(tlex)="topo"
tlex$data="lexi"
tmor<-as.data.frame(mor[which(mor == "abba" | mor == "baba" )])
colnames(tmor)="topo"
tmor$data="morpho"
tpho<-as.data.frame(pho[which(pho == "abba" | pho == "baba" )])
colnames(tpho)="topo"
tpho$data="phono"
topo<-rbind(tlex,tmor,tpho)
stopo<-topo[order(topo$data,topo$topo),]
write.table(stopo,file=paste(file,"_traits.txt", sep=""),col.names=F, row.names=T)


dlexi<-table(lex)
dmorph<-table(mor)
dphono<-table(pho)

library(gtools)
temp<-smartbind(t(as.data.frame(row.names = 1, dlexi)),
t(as.data.frame(row.names = 1, dmorph)),
t(as.data.frame(row.names = 1, dphono)),fill=0)
pat<-colnames(temp)
temp$P1<-ord[1]
temp$P2<-ord[2]
temp$P3<-ord[3]
temp$O<-ord[4]
temp$data<-c("lexi","morpho","phono")
temp<-temp[,c("P1","P2","P3","O","data",pat)]
temp$D<-(temp$abba-temp$baba)/(temp$abba+temp$baba)
temp$Dall=(sum(temp$abba)-sum(temp$baba))/(sum(temp$abba)+sum(temp$baba))
write.table(temp,file=paste(file,"_res.txt", sep=""),col.names=T)
