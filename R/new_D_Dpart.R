#!/usr/bin/env Rscript
# Théo

require('ape')
require('phangorn')


sitesquintet<-function(site){
  if (sum(site) == 1){
    if (site[1]==1) {return("baaaa")}
    else if(site[2]==1) {return("abaaa")}
    else if(site[3]==1) {return("aabaa")}
    else if(site[4]==1) {return("aaaba")}
  } else if (sum(site) == 2){
    if (site[1]+site[3] == 2){return("babaa")}
    else if (site[1]+site[4] == 2){return("baaba")}
    else if (site[2]+site[3] == 2){return("abbaa")}
    else if (site[2]+site[4] == 2){return("ababa")}
  } else if (sum(site) == 3){
    if (site[1]==0) {return("abbba")}
    else if(site[2]==0) {return("babba")}
    else if(site[3]==0) {return("bbaba")}
    else if(site[4]==0) {return("bbbaa")}
  }
}

inverser<-function(x){
  if (x[length(x)] == 1){
    return(matrix <- +(!x))
  } else {
    return(x)
  }
}

filter<-function(x){
  if (sum(x) == 1 & x[5] == 1){
    return()
  } else if (sum(x) == 4 & x[5] == 0){
    return()
  } else {return(x)}
}

Dfoiler<-function(top, site){

  names<-c("aaaba","aabaa","abaaa","baaaa","babaa","baaba","abbaa","ababa","bbbaa","bbaba","babba","abbba")
  pat<-as.data.frame(names)
  colnames(pat)<-"Var1"
  data<-as.data.frame(t(do.call("rbind",lapply(site[top,],function(x) filter(x)))))
  data<-as.data.frame(lapply(data,function(x) inverser(x)))
  p<-as.data.frame(table(unlist(lapply(data, function(x) sitesquintet(x)))))
  d<-t(merge(p,pat,all=TRUE))
  d[is.na(d)]<-0
  colnames(d)<-as.matrix(d[1, ])
  d<-as.data.frame(t(d[-1, ]))
  d<-lapply(d, function(x) as.numeric(as.character(x)))

  Dfo<-((d$babaa + d$bbbaa + d$ababa + d$aaaba) - (d$baaba + d$bbaba + d$abbaa + d$aabaa)) /
        ((d$babaa + d$bbbaa + d$ababa + d$aaaba) + (d$baaba + d$bbaba + d$abbaa + d$aabaa))
  PDfo<-binom.test(c((d$babaa + d$bbbaa + d$ababa + d$aaaba), (d$baaba + d$bbaba + d$abbaa + d$aabaa)), p = 0.5)$p.value

  Dil<-((d$abbaa + d$bbbaa + d$baaba + d$aaaba) - (d$ababa + d$bbaba + d$babaa + d$aabaa)) /
        ((d$abbaa + d$bbbaa + d$baaba + d$aaaba) + (d$ababa + d$bbaba + d$babaa + d$aabaa))
  PDil<-binom.test(c((d$abbaa + d$bbbaa + d$baaba + d$aaaba), (d$ababa + d$bbaba + d$babaa + d$aabaa)), p = 0.5)$p.value

  Dfi<-((d$babaa + d$babba + d$ababa + d$abaaa) - (d$abbaa + d$abbba + d$baaba + d$baaaa)) /
        ((d$babaa + d$babba + d$ababa + d$abaaa) + (d$abbaa + d$abbba + d$baaba + d$baaaa))
  PDfi<-binom.test(c((d$babaa + d$babba + d$ababa + d$abaaa), (d$abbaa + d$abbba + d$baaba + d$baaaa)), p = 0.5)$p.value

  Dol<-((d$baaba + d$babba + d$abbaa + d$abaaa) - (d$ababa + d$abbba + d$babaa + d$baaaa)) /
        ((d$baaba + d$babba + d$abbaa + d$abaaa) + (d$ababa + d$abbba + d$babaa + d$baaaa))
  PDol<-binom.test(c((d$baaba + d$babba + d$abbaa + d$abaaa), (d$ababa + d$abbba + d$babaa + d$baaaa)), p = 0.5)$p.value

  res<-c("P1"=top[1],
        "P2"=top[2],
        "P3"=top[3],
        "P4"=top[4],
        "P5"=top[5],
        "Dfo"=Dfo, "Pvalue_Dfo"=PDfo,
        "Dil"=Dil, "Pvalue_Dil"=PDil,
        "Dfi"=Dfi, "Pvalue_Dfi"=PDfi,
        "Dol"=Dol, "Pvalue_Dol"=PDol)

  return(res)
}

validateandreorderD<-function(arr, dist) {
  # created by Damien
  submat<-round(dist[arr,arr],5)
  if (sum(submat==max(submat))==6) {
    diag(submat)<-1
    return(names(sort(apply(submat,1,prod))))
  }
}

getquartet<-function(tr) {
  # created by Damien
  dist<-cophenetic(tr)
  allquat<-combn(tr$tip.label,4)
  RES<-do.call(rbind,apply(allquat, 2, function(x) validateandreorderD(x, dist)))
  return(RES)
}

validateandreorderDoil<-function(arr, dist) {
  submat<-dist[arr,arr]
  if (round(submat[arr[1],arr[2]],4) != round(submat[arr[3],arr[1]],4) && round(submat[arr[1],arr[3]],4) != round(submat[arr[3],arr[4]],4)){
    if (round(submat[arr[1],arr[2]],4) >= round(submat[arr[3],arr[4]],4)){
      return(c(arr[c(3,4)],arr[c(1,2)],arr[5]))
    } else {
      return(arr)
    }
  }
}

getquintet<-function(top){
  dist<-cophenetic(tr)
  duos<-unique(top[,c(1,2)])
  tres<-unique(top[,c(1,2,4)])
  mylist <- list()
  for (j in 1:nrow(tres)){
    mylist[[j]]<-data.frame(matrix(NA, ncol = 5))
    a=tres[j,]
    list<-top[top[,1]%in%a[1] & top[,2]%in%a[2] & top[,4]%in%a[3],3]
    second_pair<-top[top[,1]%in%list & top[,2]%in%list & top[,3]%in%a & top[,4]%in%a[3],]
    pair<-unique(second_pair[,c(1,2)])
    if (nrow(pair) != 0){
      for (i in 1:nrow(pair)){
        mylist[[j]][i,]<-(c(a[1],a[2],pair[i,1],pair[i,2],a[3]))
      }
    }
  }
  df<-unique(na.omit(do.call("rbind",mylist)))
  res<-do.call("rbind",apply(df,1,function(x) validateandreorderDoil(unlist(x), dist)))
  return(unique(res))
}


sound<-read.table('data3', h = T)
tr<-read.tree('spe_tree')

spnd<-c(tr$tip.label, tr$node.label)

pretopologiesD<-getquartet(tree)

top=pretopologiesD


(pretopologiesDfoil<-getquintet(pretopologiesD))


(a=Dfoiler(topologiesDfoil[1,],sound))
