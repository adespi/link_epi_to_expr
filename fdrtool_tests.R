setwd("~/mit_internship/link_epi_to_expr")
library("fdrtool")
file<-'correlations_para/correlations_17_17740325_SREBF1.csv.gz'
#file<-'correlations_para/correlations_2_127864931_BIN1.csv.gz'
csv<-read.csv(file,head=T)
file<-gsub("[/.]", "_", file)
gene<-strsplit(file, '_')[[1]][6]
row.names(csv)<-paste(cbind(csv[1],seq(919))$X,cbind(csv[1],seq(919))$`seq(919)`)
csv[1]<- list(NULL)
csv[is.na(csv)] <- 0

p<-0.01

fdr.out = fdrtool(as.vector(t(csv[1,])), statistic = "correlation")
plot(fdr.out$pval)#,ylim=c(0, 0.05))
title(main = paste(gene,"p-values for one mark"))
plot(as.numeric(csv[1,]))
a<-which(fdr.out$pval<p,TRUE)
b<-as.numeric(csv[1,])[fdr.out$pval<p]
points(a,b,pch=19,col='red')
title(main = paste(gene," corr coef for one mark (",length(a)," points with p<",p,")",sep=""))

'csv2<-apply((csv),2,max)
plot(csv2)
fdr.out = fdrtool(csv2, statistic = "correlation")
plot(fdr.out$pval)'

csv2<-apply(abs(csv),2,max)
fdr.out = fdrtool(csv2, statistic = "correlation")
plot(fdr.out$pval)
title(main = paste(gene,"p-values for max(mark)"))
plot(csv2)
a<-which(fdr.out$pval<p,TRUE)
b<-csv2[fdr.out$pval<p]
points(a,b,pch=19,col='red')
title(main = paste(gene," corr coef for max(mark) (",length(a)," points with p<",p,")",sep=""))


'#to keep max(csv) but with negative values too
a<-apply(abs(csv),2,which.max)
x<-1
for (i in a){
  a[x]<-csv[i,x]
  x<-x+1
}
a
plot(a)
fdr.out = fdrtool(a, statistic = "correlation")
plot(fdr.out$pval)'


#plot(fdr.out$pval<0.00001)