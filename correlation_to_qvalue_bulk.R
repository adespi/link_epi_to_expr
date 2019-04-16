library("fdrtool")
p <- 0.05
#file<-'correlations/correlations_17_17740325_SREBF1.csv.gz'
#file<-'correlations/correlations_2_127864931_BIN1.csv.gz'
for (file in Sys.glob(paste("correlations/*.csv.gz", sep = ''))) {
  print(file)
  csv <- read.csv(file, head = T, check.names = F)
  gene <- tail(strsplit(gsub("[/.]", "_", file), "_")[[1]],n=3)[1]
  row.names(csv) <-
    paste(cbind(csv[1], seq(919))$X, cbind(csv[1], seq(919))$`seq(919)`)
  csv[1] <- list(NULL)
  csv[is.na(csv)] <- 0
  csv[apply(abs(csv), 2, max)==0] <- list(NULL)
  
  fdr.out = fdrtool(
    as.vector(as.matrix(csv)),
    statistic = "correlation",
    plot = F,
    verbose = F
  )
  
  #min(abs(csv)[matrix(fdr.out$qval,919,5000)<0.05])
  qvals = data.frame(matrix(fdr.out$qval,919,length(colnames(csv))))
  rownames(qvals) = rownames(csv)
  colnames(qvals) = colnames(csv)
  qvals = apply(qvals, 2, min)
  qvals = qvals[qvals < 0.05]
  write.csv(t(qvals[qvals < 0.05]), file = gzfile(gsub("correlations/","correlations_small/", file)), row.names = FALSE)
}
