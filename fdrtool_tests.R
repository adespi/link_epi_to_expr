setwd("~/mit_internship/link_epi_to_expr")
library("fdrtool")
file<-'correlations_para/correlations_17_17740325_SREBF1.csv.gz'
file<-'correlations_para/correlations_2_127864931_BIN1.csv.gz'
csv<-read.csv(file,head=T)
file<-gsub("[/.]", "_", file)
gene<-strsplit(file, '_')[[1]][6]
row.names(csv)<-paste(cbind(csv[1],seq(919))$X,cbind(csv[1],seq(919))$`seq(919)`)
csv[1]<- list(NULL)
csv[is.na(csv)] <- 0

p<-0.05
qvals=data.frame(matrix(ncol=5000,nrow=919))
for (x in 1:919) {
  fdr.out = fdrtool(as.vector(t(csv[x,])), statistic = "correlation",plot = F,verbose = F)
  qvals[x,]=fdr.out$qval
  print(x)
  #plot(fdr.out$qval<0.05)#,ylim=c(0, 0.05))
  #title(main = paste(gene,"q-values for one mark :",x))
}
#qvals=qvals[sample(nrow(qvals)),]
coefs=c(1)
corr=cor(t(qvals))
sd(qvals[1,])
min(qvals)
'#corr[is.na(corr)]=0
corr=as.vector(corr)
corr=corr[corr!=1]
plot(sort(corr))
median(corr)
mean(corr)
min(corr)
max(corr)'
for (x in 2:919) {
  f=1
  for (y in (1:(x-1))){
    f=f- abs(corr[y, x])*coefs[y]
  }
  coefs=append(coefs,max(0,f))
  print(x)
}
plot(coefs)
library("matrixStats", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
qvl=as.matrix((qvals^(coefs)))
qvl[is.na(qvl)]=1
qvl=colProds(qvl)
#qvl[1:15]
plot(qvl)
#qvl1=qvl
#k=rownames(as.data.frame(sort(apply(qvals, 2,mean))))
k=as.integer(gsub("X","",rownames(as.data.frame(sort(apply(qvals, 2,mean))))))
plot(qvl,qvl1)
plot(qvl,as.matrix(apply(qvals, 2,min)))
###plot(as.matrix(apply(qvals, 2,min)),as.matrix(apply(qvals<p, 2,sum)),xlim=c(0,0.05))
plot(as.matrix(apply(qvals, 2,min)),as.matrix(apply(qvals<p, 2,sum)),xlim=c(0,0.05))
good_points=as.matrix(apply(qvals<p, 2,sum))>0
sum(good_points)
points(as.matrix(apply(qvals, 2,min))[good_points],as.matrix(apply(qvals<p, 2,sum))[good_points],col='red',pch=19)
plot(qvl,log(as.matrix(apply(qvals<p, 2,sum))))
plot(apply(qvals, 2,mean))
plot(apply(qvals, 2,min))
plot(apply(qvals, 1,min))
plot(apply(qvals<p, 2,sum))
abline(h=sum(qvl<p))


plot()

#colProds(matrix(1:12, nrow = 3, ncol = 4)^c(1,2,3))

""


write.csv(qvals, file = "qvals.csv")
plot(apply(qvals,1,min))
which.max(apply(qvals<0.05, 2,sum))
plot(apply(qvals<0.05, 1,sum))
cols2 <- colorRampPalette(c("blue","white"))(256)
image(as.matrix(t(qvals)),  col = cols2, zlim=c(0, 0.05))
title(main="image of the qvalues",
      xlab="posiyions", ylab="histone marks") 
sum(qvals[714,]<0.05)
which.max(apply(qvals<0.05, 0,sum))

fdr.out = fdrtool(as.vector(t(csv[3,])), statistic = "correlation",plot = F)
'plot(fdr.out$pval)#,ylim=c(0, 0.05))
title(main = paste(gene,"p-values for one mark"))'
plot(fdr.out$qval)#,ylim=c(0, 0.05))
title(main = paste(gene,"q-values for one mark"))
'plot(as.numeric(csv[1,]))
a<-which(fdr.out$qval<p,TRUE)
b<-as.numeric(csv[1,])[fdr.out$qval<p]
points(a,b,pch=19,col="red")
title(main = paste(gene," corr coef for one mark (",length(a)," points with p<",p,")",sep=""))
'
'csv2<-apply((csv),2,max)
plot(csv2)
fdr.out = fdrtool(csv2, statistic = "correlation")
plot(fdr.out$pval)

csv2<-apply(abs(csv),2,max)
fdr.out = fdrtool(csv2, statistic = "correlation")
plot(fdr.out$pval)
title(main = paste(gene,"p-values for max(mark)"))
plot(fdr.out$qval)
title(main = paste(gene,"q-values for max(mark)"))
plot(csv2)
a<-which(fdr.out$qval<p,TRUE)
b<-csv2[fdr.out$qval<p]
points(a,b,pch=19,col="red")
title(main = paste(gene," corr coef for max(mark) (",length(a)," points with p<",p,")",sep=""))



#to keep max(csv) but with negative values too
a<-apply(abs(csv),2,which.max)
x<-1
for (i in a){
  a[x]<-csv[i,x]
  x<-x+1
}
a
plot(a)
fdr.out = fdrtool(a, statistic = "correlation")
plot(fdr.out$qval)


#plot(fdr.out$pval<0.00001)'