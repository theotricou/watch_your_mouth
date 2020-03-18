


require('ape')
require('phangorn')

tree = read.tree('spe_gra')
# plot(tree, show.node = T)

tr = read.table( "signi_ancestor_internal_ABBA_BABA.txt", sep = "\t", h = T)

par(new=TRUE)


tttt = read.tree('spe_dist')
ttttspnd<-c(tttt$tip.label, tttt$node.label)
spnd<-c(tree$tip.label, tree$node.label)

plot(tree,
  no.margin = TRUE,
  edge.width = 2,
  direction = "downwards",
  type = "phylogram",
)
nodelabels(spnd[c(29:45)], adj = c(0.5, -0.8), cex = 0.5)


trcol <- function(col, alpha=1){
    if(missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2,
          function(x)
              rgb(x[1], x[2], x[3], alpha=alpha))}




edge<-tree$edge
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
x<-lastPP$xx
y<-lastPP$yy
par(new=TRUE)


tr = read.table( "signi_ancestor_internal_ABBA_BABA.txt", sep = "\t", h = T)
tr$P1 = as.character(tr$P1)
tr$P2 = as.character(tr$P2)
for (i in 1:nrow(tr)) {
  if (as.numeric(as.character(tr[i,10])) > 0) {
    temp = tr[i,1]
    tr[i,1] = tr[i,2]
    tr[i,2] = temp
  }
}


tr = tr[order(tr$P2),]


n11 = tr[tr$P2 == "N11",]
n15 = tr[tr$P2 == "N15",]
n9 = tr[tr$P2 == "N9",]

pdf("glotto_tree_P3tr_pdf",width=29.7/2.54, height=21/2.54)

par(new=TRUE)

plot(tree,
  show.node = T,
  no.margin = TRUE,
  edge.width = 2,
  direction = "downwards",
  type = "phylogram",
  srt = 0,
  adj = -0.8
)
title("Tranfers events from P3 to P1/P2", line = -2)


# P3 to P12

matP3 <- matrix(ncol=5, nrow=nrow(tr))
s <- seq(nrow(tr))  # one shorter than data
for (i in 1:nrow(tr)) {
  if (tr[i,10] < 0) {
    from<-tr[i,3]
    to<-tr[i,1]
  } else {
    from<-tr[i,3]
    to<-tr[i,2]
  }
  wherefrom<-which(spnd==from)
  whereto<-which(spnd==to)
  fromdad<-edge[edge[,2]==wherefrom,1]
  todad<-edge[edge[,2]==whereto,1]
  # ywherefrom <- y[wherefrom] - (y[wherefrom] - y[fromdad])/2
  ywhereto <- y[whereto] - (y[whereto] - y[todad])/2
  ywherefrom <- rnorm(1, mean= y[wherefrom] - (y[wherefrom] - y[fromdad])/2)

  matP3[i,] <- c(x[wherefrom], ywherefrom,x[whereto],ywhereto,1)
}
arrows(as.numeric(matP3[s,1]),as.numeric(matP3[s,2]),as.numeric(matP3[s,3]),as.numeric(matP3[s,4]), col ="red",lwd=3, length = 0.2)


dev.off()



pdf("glotto_tree_P4tr_pdf",width=29.7/2.54, height=21/2.54)

par(new=TRUE)

plot(tree,
  show.node = T,
  no.margin = TRUE,
  edge.width = 2,
  direction = "downwards",
  type = "phylogram",
  srt = 0,
  adj = -0.8
)
title("Tranfers events from P4 to P1/P2", line = -2)



# P4 to P12
matP4 <- matrix(ncol=5, nrow=nrow(tr))
s <- seq(nrow(tr))  # one shorter than data
for (i in 1:nrow(tr)) {
  if (tr[i,10] < 0) {
    from<-tr[i,4]
    to<-tr[i,2]
  } else {
    from<-tr[i,4]
    to<-tr[i,1]
  }
  wherefrom<-which(spnd==from)
  whereto<-which(spnd==to)
  fromdad<-edge[edge[,2]==wherefrom,1]
  todad<-edge[edge[,2]==whereto,1]
  # ywherefrom <- y[wherefrom] - (y[wherefrom] - y[fromdad])/2
  ywhereto <- y[whereto] - (y[whereto] - y[todad])/2

  # ywherefrom <- runif(1, min = y[wherefrom],  max = y[wherefrom] - (y[wherefrom] - y[fromdad]))
  ywherefrom <- rnorm(1, (mean= y[wherefrom] - (y[wherefrom] - y[fromdad])/2))

  matP4[i,] <- c(x[wherefrom], ywherefrom,x[whereto],ywhereto,1)
}

arrows(as.numeric(matP4[s,1]),as.numeric(matP4[s,2]),as.numeric(matP4[s,3]),as.numeric(matP4[s,4]), col ="blue",lwd=3, length = 0.2)

dev.off()








tr = read.table( "signi_ancestor_internal_ABBA_BABA.txt", sep = "\t", h = T)
tr$P1 = as.character(tr$P1)
tr$P2 = as.character(tr$P2)
for (i in 1:nrow(tr)) {
  if (as.numeric(as.character(tr[i,10])) > 0) {
    print(i)
    temp = tr[i,1]
    tr[i,1] = tr[i,2]
    tr[i,2] = temp
  }
}
tr = tr[order(tr$P2),]
n11 = tr[tr$P2 == "N11",]
n15 = tr[tr$P2 == "N15",]
n9 = tr[tr$P2 == "N9",]
nmara = tr[tr$P2 == "mara1378",]
nsout = tr[tr$P2 == "sout2672",]

pdf("glotto_N11.pdf",width=29.7/2.54, height=21/2.54)
tr = n11
par(new=TRUE)
plot(tree,
  no.margin = TRUE,
  edge.width = 2,
  direction = "downwards",
  type = "phylogram",
)
title("ABBA BABA tranfers events relative to P2 = N11", line = -2, adj = 0.1)
matP4 <- matrix(ncol=5, nrow=nrow(tr))
s <- seq(nrow(tr))  # one shorter than data
for (i in 1:nrow(tr)) {
  from<-tr[i,4]
  to<-tr[i,2]
  wherefrom<-which(spnd==from)
  whereto<-which(spnd==to)
  fromdad<-edge[edge[,2]==wherefrom,1]
  todad<-edge[edge[,2]==whereto,1]
  ywherefrom <- rnorm(1, mean= y[wherefrom] - (y[wherefrom] - y[fromdad])/2, 2)
  ywhereto <- rnorm(1, mean= y[whereto] - (y[whereto] - y[todad])/2, 2)
  matP4[i,] <- c(x[wherefrom], ywherefrom,x[whereto],ywhereto,1)}
matP3 <- matrix(ncol=5, nrow=nrow(tr))
s <- seq(nrow(tr))  # one shorter than data
for (i in 1:nrow(tr)) {
  from<-tr[i,3]
  to<-tr[i,1]
  wherefrom<-which(spnd==from)
  whereto<-which(spnd==to)
  fromdad<-edge[edge[,2]==wherefrom,1]
  todad<-edge[edge[,2]==whereto,1]
  ywhereto <- y[whereto] - (y[whereto] - y[todad])/2
  ywherefrom <- rnorm(1, mean= y[wherefrom] - (y[wherefrom] - y[fromdad])/2, 2)
  ywhereto <- rnorm(1, mean= y[whereto] - (y[whereto] - y[todad])/2, 2)
  matP3[i,] <- c(x[wherefrom], ywherefrom,x[whereto],ywhereto,1)}
arrows(as.numeric(matP4[s,1]),as.numeric(matP4[s,2]),as.numeric(matP4[s,3]),as.numeric(matP4[s,4]), col ="red",lwd=1.5, length = 0.15, code = 3)
arrows(as.numeric(matP3[s,1]),as.numeric(matP3[s,2]),as.numeric(matP3[s,3]),as.numeric(matP3[s,4]), col ="blue",lwd=1.5, length = 0.15, code = 3)
nodelabels(ttttspnd[c(29:45)], adj = c(0.5, -0.8), cex = 0.7)

dev.off()




pdf("languagetree.pdf",width=29.7/2.54, height=21/2.54)


plot(tree,
  no.margin = TRUE,
  edge.width = 2,
  direction = "downwards",
  type = "phylogram",
)
title("Arbre des langues indo-aryennes", line = -2, adj = 0.1)
nodelabels(ttttspnd[c(29:45)], adj = c(0.5, -0.8), cex = 0.7)



dev.off()



pdf("glotto_N15.pdf",width=29.7/2.54, height=21/2.54)
tr = n15
par(new=TRUE)
plot(tree,
  no.margin = TRUE,
  edge.width = 2,
  direction = "downwards",
  type = "phylogram",
)
title("ABBA BABA tranfers events relative to P2 = Mewati-Gojri", line = -2, adj = 0.1)
matP4 <- matrix(ncol=5, nrow=nrow(tr))
s <- seq(nrow(tr))  # one shorter than data
for (i in 1:nrow(tr)) {
  from<-tr[i,4]
  to<-tr[i,2]
  wherefrom<-which(spnd==from)
  whereto<-which(spnd==to)
  fromdad<-edge[edge[,2]==wherefrom,1]
  todad<-edge[edge[,2]==whereto,1]
  ywherefrom <- rnorm(1, mean= y[wherefrom] - (y[wherefrom] - y[fromdad])/2, 2)
  ywhereto <- rnorm(1, mean= y[whereto] - (y[whereto] - y[todad])/2, 2)
  matP4[i,] <- c(x[wherefrom], ywherefrom,x[whereto],ywhereto,1)}
matP3 <- matrix(ncol=5, nrow=nrow(tr))
s <- seq(nrow(tr))  # one shorter than data
for (i in 1:nrow(tr)) {
  from<-tr[i,3]
  to<-tr[i,1]
  wherefrom<-which(spnd==from)
  whereto<-which(spnd==to)
  fromdad<-edge[edge[,2]==wherefrom,1]
  todad<-edge[edge[,2]==whereto,1]
  ywhereto <- y[whereto] - (y[whereto] - y[todad])/2
  ywherefrom <- rnorm(1, mean= y[wherefrom] - (y[wherefrom] - y[fromdad])/2, 2)
  ywhereto <- rnorm(1, mean= y[whereto] - (y[whereto] - y[todad])/2, 2)
  matP3[i,] <- c(x[wherefrom], ywherefrom,x[whereto],ywhereto,1)}
arrows(as.numeric(matP4[s,1]),as.numeric(matP4[s,2]),as.numeric(matP4[s,3]),as.numeric(matP4[s,4]), col ="red",lwd=1.5, length = 0.15, code = 3)
arrows(as.numeric(matP3[s,1]),as.numeric(matP3[s,2]),as.numeric(matP3[s,3]),as.numeric(matP3[s,4]), col ="blue",lwd=1.5, length = 0.15, code = 3)
nodelabels(ttttspnd[c(29:45)], adj = c(0.5, -0.8), cex = 0.7)

dev.off()



pdf("glotto_N9.pdf",width=29.7/2.54, height=21/2.54)
tr = n9
par(new=TRUE)
plot(tree,
  no.margin = TRUE,
  edge.width = 2,
  direction = "downwards",
  type = "phylogram",
)
title("ABBA BABA tranfers events relative to P2 = Maithili-Magahi", line = -2, adj = 0.1)
matP4 <- matrix(ncol=5, nrow=nrow(tr))
s <- seq(nrow(tr))  # one shorter than data
for (i in 1:nrow(tr)) {
  from<-tr[i,4]
  to<-tr[i,2]
  wherefrom<-which(spnd==from)
  whereto<-which(spnd==to)
  fromdad<-edge[edge[,2]==wherefrom,1]
  todad<-edge[edge[,2]==whereto,1]
  ywherefrom <- rnorm(1, mean= y[wherefrom] - (y[wherefrom] - y[fromdad])/2, 2)
  ywhereto <- rnorm(1, mean= y[whereto] - (y[whereto] - y[todad])/2, 2)
  matP4[i,] <- c(x[wherefrom], ywherefrom,x[whereto],ywhereto,1)}
matP3 <- matrix(ncol=5, nrow=nrow(tr))
s <- seq(nrow(tr))  # one shorter than data
for (i in 1:nrow(tr)) {
  from<-tr[i,3]
  to<-tr[i,1]
  wherefrom<-which(spnd==from)
  whereto<-which(spnd==to)
  fromdad<-edge[edge[,2]==wherefrom,1]
  todad<-edge[edge[,2]==whereto,1]
  ywhereto <- y[whereto] - (y[whereto] - y[todad])/2
  ywherefrom <- rnorm(1, mean= y[wherefrom] - (y[wherefrom] - y[fromdad])/2, 2)
  ywhereto <- rnorm(1, mean= y[whereto] - (y[whereto] - y[todad])/2, 2)
  matP3[i,] <- c(x[wherefrom], ywherefrom,x[whereto],ywhereto,1)}
arrows(as.numeric(matP4[s,1]),as.numeric(matP4[s,2]),as.numeric(matP4[s,3]),as.numeric(matP4[s,4]), col ="red",lwd=1.5, length = 0.15, code = 3)
arrows(as.numeric(matP3[s,1]),as.numeric(matP3[s,2]),as.numeric(matP3[s,3]),as.numeric(matP3[s,4]), col ="blue",lwd=1.5, length = 0.15, code = 3)
nodelabels(ttttspnd[c(29:45)], adj = c(0.5, -0.8), cex = 0.7)

dev.off()






pdf("glotto_mara.pdf",width=29.7/2.54, height=21/2.54)
tr = nmara
par(new=TRUE)
plot(tree,show.node = T, no.margin = TRUE, edge.width = 2, direction = "downwards",
  type = "phylogram", srt = 0, adj = -0.8)
title("ABBA BABA tranfers events relative to P2 = mara", line = -2)
matP4 <- matrix(ncol=5, nrow=nrow(tr))
s <- seq(nrow(tr))  # one shorter than data
for (i in 1:nrow(tr)) {
  from<-tr[i,4]
  to<-tr[i,2]
  wherefrom<-which(spnd==from)
  whereto<-which(spnd==to)
  fromdad<-edge[edge[,2]==wherefrom,1]
  todad<-edge[edge[,2]==whereto,1]
  ywherefrom <- rnorm(1, mean= y[wherefrom] - (y[wherefrom] - y[fromdad])/2, 2)
  ywhereto <- rnorm(1, mean= y[whereto] - (y[whereto] - y[todad])/2, 2)
  matP4[i,] <- c(x[wherefrom], ywherefrom,x[whereto],ywhereto,1)}
matP3 <- matrix(ncol=5, nrow=nrow(tr))
s <- seq(nrow(tr))  # one shorter than data
for (i in 1:nrow(tr)) {
  from<-tr[i,3]
  to<-tr[i,1]
  wherefrom<-which(spnd==from)
  whereto<-which(spnd==to)
  fromdad<-edge[edge[,2]==wherefrom,1]
  todad<-edge[edge[,2]==whereto,1]
  ywhereto <- y[whereto] - (y[whereto] - y[todad])/2
  ywherefrom <- rnorm(1, mean= y[wherefrom] - (y[wherefrom] - y[fromdad])/2, 2)
  ywhereto <- rnorm(1, mean= y[whereto] - (y[whereto] - y[todad])/2, 2)
  matP3[i,] <- c(x[wherefrom], ywherefrom,x[whereto],ywhereto,1)}
arrows(as.numeric(matP4[s,1]),as.numeric(matP4[s,2]),as.numeric(matP4[s,3]),as.numeric(matP4[s,4]), col ="red",lwd=1.5, length = 0.15)
arrows(as.numeric(matP3[s,1]),as.numeric(matP3[s,2]),as.numeric(matP3[s,3]),as.numeric(matP3[s,4]), col ="blue",lwd=1.5, length = 0.15)
dev.off()
pdf("glotto_sout.pdf",width=29.7/2.54, height=21/2.54)
tr = nsout
par(new=TRUE)
plot(tree,show.node = T, no.margin = TRUE, edge.width = 2, direction = "downwards",
  type = "phylogram", srt = 0, adj = -0.8)
title("ABBA BABA tranfers events relative to P2 = sout", line = -2)
matP4 <- matrix(ncol=5, nrow=nrow(tr))
s <- seq(nrow(tr))  # one shorter than data
for (i in 1:nrow(tr)) {
  from<-tr[i,4]
  to<-tr[i,2]
  wherefrom<-which(spnd==from)
  whereto<-which(spnd==to)
  fromdad<-edge[edge[,2]==wherefrom,1]
  todad<-edge[edge[,2]==whereto,1]
  ywherefrom <- rnorm(1, mean= y[wherefrom] - (y[wherefrom] - y[fromdad])/2, 2)
  ywhereto <- rnorm(1, mean= y[whereto] - (y[whereto] - y[todad])/2, 2)
  matP4[i,] <- c(x[wherefrom], ywherefrom,x[whereto],ywhereto,1)}
matP3 <- matrix(ncol=5, nrow=nrow(tr))
s <- seq(nrow(tr))  # one shorter than data
for (i in 1:nrow(tr)) {
  from<-tr[i,3]
  to<-tr[i,1]
  wherefrom<-which(spnd==from)
  whereto<-which(spnd==to)
  fromdad<-edge[edge[,2]==wherefrom,1]
  todad<-edge[edge[,2]==whereto,1]
  ywhereto <- y[whereto] - (y[whereto] - y[todad])/2
  ywherefrom <- rnorm(1, mean= y[wherefrom] - (y[wherefrom] - y[fromdad])/2, 2)
  ywhereto <- rnorm(1, mean= y[whereto] - (y[whereto] - y[todad])/2, 2)
  matP3[i,] <- c(x[wherefrom], ywherefrom,x[whereto],ywhereto,1)}
arrows(as.numeric(matP4[s,1]),as.numeric(matP4[s,2]),as.numeric(matP4[s,3]),as.numeric(matP4[s,4]), col ="red",lwd=1.5, length = 0.15)
arrows(as.numeric(matP3[s,1]),as.numeric(matP3[s,2]),as.numeric(matP3[s,3]),as.numeric(matP3[s,4]), col ="blue",lwd=1.5, length = 0.15)
dev.off()
