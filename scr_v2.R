
d = read.table('cognate_matrix', h=T)
patterner <- function(x){
  if (length(table(x)) == 2){
    if(table(x)[1] == 2){
      return(x)
    }
  }
}


a = na.omit(d[,c("deu","eng","fra","bul")])
res = do.call('rbind', apply(a, 1, function(x) patterner(x)))

abba=0
baba=0
for (i in 1:nrow(res)){
  if (res[i,1] == res[i,4]){
    abba = abba + 1
  }else if (res[i,1] == res[i,3]){
    baba = baba + 1
  }
}

abba
baba

a = na.omit(d[,c("ron","ita","bul", "asm")])
