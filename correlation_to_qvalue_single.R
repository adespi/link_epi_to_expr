args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please specify a file", call. = FALSE)
}
#file<-'correlations_para/correlations_17_17740325_SREBF1.csv.gz'
#file<-'correlations_para/correlations_2_127864931_BIN1.csv.gz'
file <- args[1]
library("fdrtool")
p <- 0.05

print(file)
csv <- read.csv(file, head = T, check.names = F)
gene <- strsplit(gsub("[/.]", "_", file), "_")[[1]][6]
row.names(csv) <-
  paste(cbind(csv[1], seq(919))$X, cbind(csv[1], seq(919))$`seq(919)`)
csv[1] <- list(NULL)
csv[is.na(csv)] <- 0

qvals = data.frame(matrix(ncol = length(colnames(csv)), nrow = 919))
rownames(qvals) = rownames(csv)
colnames(qvals) = colnames(csv)
for (x in 1:919) {
  fdr.out = fdrtool(
    as.vector(t(csv[x, ])),
    statistic = "correlation",
    plot = F,
    verbose = F
  )
  qvals[x, ] = fdr.out$qval
  #plot(fdr.out$qval<0.05)#,ylim=c(0, 0.05))
  #title(main = paste(gene,"q-values for one mark :",x))
}
qvals = apply(qvals, 2, min)
qvals = qvals[qvals < 0.05]
write.csv(t(qvals), file = gzfile(paste(
  "correlations_small/", basename(file), sep = ""
)), row.names = FALSE)

warnings()