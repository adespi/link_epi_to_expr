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

plot(as.numeric(csv[1,]))
title(main = paste(gene,"corr coef for one mark"))
fdr.out = fdrtool(as.vector(t(csv[1,])), statistic = "correlation")
plot(fdr.out$pval)#,ylim=c(0, 0.05))
title(main = paste(gene,"p-values for one mark"))


'a<-apply((csv),2,max)
plot(a)
fdr.out = fdrtool(a, statistic = "correlation")
plot(fdr.out$pval)'

a<-apply(abs(csv),2,max)
plot(a)
title(main = paste(gene,"corr coef for max(mark)"))
fdr.out = fdrtool(a, statistic = "correlation")
plot(fdr.out$pval)
title(main = paste(gene,"p-values for max(mark)"))

a<-apply(abs(csv),2,which.max)
x<-1
for (i in a){
  a[x]<-csv[i,x]
  x<-x+1
}
a
plot(a)
fdr.out = fdrtool(a, statistic = "correlation")
plot(fdr.out$pval)


plot(fdr.out$pval<0.00001)