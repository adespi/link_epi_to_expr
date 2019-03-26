#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Please specify a gene number", call.=FALSE)
}
#setwd("~/link_epi_to_expr")
#beware suppressWarnings applied to correlation
chromosome<-args[1]
gene<-args[2]
gene_name<-args[3]
#gene<-strsplit(gsub("[a-z,/,.]","\\1",Sys.glob("temp/predictions_with_names/*"))[1],'_')[[1]][3]
i=1
list_positions<-c()
for (file in gsub("[a-z,/,.]","\\1",Sys.glob(paste("temp/",chromosome,'_',gene,"/predictions_with_names/predictions",chromosome,'_',gene,"_*",sep='')))){
    list_positions[i]<-strsplit(gsub("[a-z,/,.]","\\1",Sys.glob(paste("temp/",chromosome,'_',gene,"/predictions_with_names/predictions",chromosome,'_',gene,"_*",sep=''))),'_')[[i]][6]
    i<-i+1
}
list_positions<-as.character(format(as.numeric(sort(as.numeric(list_positions))), scientific=F,trim=T))
expression <- read.table(paste("temp/",chromosome,'_',gene,"/expression/",chromosome,'_',gene,".tsv",sep=''), head=T)

for (position in list_positions){
    predictions <- read.table(paste("temp/",chromosome,'_',gene,"/predictions_with_names/predictions",chromosome,'_',gene,"_",position,".fa.tsv",sep=''), head=T,row.names=1)
    #create correlation file
    if (position==list_positions[1])
        correlations <- data.frame(row.names=colnames(predictions)[4:dim(predictions)[2]])
    correlations<-cbind(correlations,data.frame(g=rep(0,919)))
    colnames(correlations)[colnames(correlations)=='g']<-paste("p",position,sep='')
    #add a column for expression
    predictions$expression<-predictions$metadata.ranges.id
    a=1
    for (i in row.names(predictions)){
        predictions$expression[a]<-expression[1,match(c(substring(i,5,11)),colnames(expression))]
        a<-a+1
    }

    for (i in colnames(predictions[4:(dim(predictions)[2]-1)])){
        correlations[i,paste("p",position,sep='')]<-signif(suppressWarnings(cor(eval(parse(text = paste("predictions$",i,sep=""))),predictions$expression)),digits=3)
        #print(paste(i,cor(eval(parse(text = paste("predictions$",i,sep=""))),predictions$expression)))
        #print()
    }
}
#warnings()
#correlations<-read.csv('correlations_before_chr_change/correlations_1_53704282.csv',head=T,row.names=1,check.names=FALSE)
feature_name<-read.table('deepsea_postprocessing/feature_name',row.names=1)
for (pred in rownames(correlations)){
    #print(feature_name[row.names(feature_name)==substring(pred,7,9),])
    #break
    rownames(correlations)[rownames(correlations)==pred]<- paste(as.character(feature_name[row.names(feature_name)==substring(pred,7,9),]),"(",substring(pred,7,9),")",sep='')
}
write.csv(correlations,file=paste("correlations_para/correlations_",chromosome,'_',gene,'_',gene_name,".csv",sep=''))
#class(predictions)
#dim(predictions)
#summary(predictions)






#setwd("~/link_epi_to_expr")
#predictions <- read.table("predictions.tsv", head=T,row.names=1)
#expression <- read.table("expression/569076.tsv", head=T)
#create correlation file
#correlations <- data.frame(row.names=colnames(predictions)[5:dim(predictions)[2]])
#correlations<-cbind(correlations,data.frame(g569076=rep(0,919)))
#add a column for expression
#predictions$expression<-predictions$metadata.ranges.id
#a=1
#for (i in row.names(predictions)){
#    predictions$expression[a]<-expression[1,match(c(substring(i,4,10)),colnames(expression))]
#    a<-a+1
#}

#for (i in colnames(predictions[5:(dim(predictions)[2]-1)])){
#    correlations[i,"g569076"]<-cor(eval(parse(text = paste("predictions$",i,sep=""))),predictions$expression)
    #print(paste(i,cor(eval(parse(text = paste("predictions$",i,sep=""))),predictions$expression)))
    #print()
#}
#class(predictions)
#dim(predictions)
#summary(predictions)
