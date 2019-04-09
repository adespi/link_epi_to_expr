csv<-read.csv('correlations/correlations_5_158526788_EBF1.csv',head=T,row.names=1)
plot(apply(abs(csv[,2000:3000]),2,max))
plot(apply(abs(csv),2,max))
cols2 <- colorRampPalette(c("blue","white","red"))(256)
image(as.matrix(csv),  col = cols2, zlim=c(-0.25, 0.25))
sort(rowMeans(csv,na.rm=T))
plot(seq(1,4987),csv["K562|E2F4|None(620)",])
plot(seq(1,4987),csv["Monocytes-CD14+RO01746|H3K27ac|None(856)",])

plot(seq(1,4987),csv["Monocytes-CD14+RO01746|H3K27ac|None(856)",])
panel.smooth(seq(1,4987),as.numeric(csv["Monocytes-CD14+RO01746|H3K27ac|None(856)",]))
abline(h=0,col='yellow')

image(as.matrix(csv[c("K562|E2F4|None(620)","K562|RFX5|None(656)","K562|RBBP5|None(159)","A549|USF-1|EtOH_0.02pct(195)","K562|CTCF|None(619)"),]),  col = cols2, zlim=c(-0.25, 0.25))
order(rowMeans(csv,na.rm=T))


csv1<-read.csv('correlations/correlations_5_88199922_test_log.csv',head=T,row.names=1)
plot(apply(abs(csv1),2,max))
csv2<-read.csv('correlations/logcorrelations_5_88199922_test_log.csv',head=T,row.names=1)
plot(apply(abs(csv2),2,max))

for (position in seq(1,length(csvbin))){
   text(position, apply(abs(csvbin),2,max)[position],colnames(csvbin)[position])
}

csvbin[is.na(csvbin)] <- 0

for (file in Sys.glob(paste("correlations/*.csv.gz",sep=''))){
   csv<-read.csv(file,head=T)
   csv[1]<- list(NULL)
   csv[is.na(csv)] <- 0
   jpeg(paste(file,'.jpg',sep=''))
   plot(apply(abs(csv),2,max))
   title(basename(file))
   dev.off()
}

jpeg("correlations_graph.jpg")
#def.par<-par(no.readonly=T)
#layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),3,5,byrow=T))
par(mfrow=c(3,5))
for (file in Sys.glob(paste("correlations/*.csv",sep=''))){
   csv<-read.csv(file,head=T,row.names=1)
   plot(apply(abs(csv),2,max))
   title(basename(file))
}
dev.off()

D<-matrix(c(7,1,4,5,2,2,4,5,5,2,5,1,7,8,0,7),nrow=4,ncol=4,byrow=TRUE)
E<-D
D
cor(D)
D
for (i in seq(1,dim(D)[2])){
   D[,i]<- D[,i]-mean(D[,i])
   D[,i]<- D[,i]/norm(as.matrix(D[,i]),"f")
}
D
cor(E)
t(D)%*%D

D<-matrix(c(7,1,4,5,2,2,4,5,5,2,5,1,7,8,0,7),nrow=4,ncol=4,byrow=TRUE)
for (i in seq(100000)){
   D<-cor(E[,1],E[,2])
}

D<-matrix(c(7,1,4,5,2,2,4,5,5,2,5,1,7,8,0,7),nrow=4,ncol=4,byrow=TRUE)
for (i in seq(100000)){
   D<-t(D)%*%D
}

