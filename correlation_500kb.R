setwd("~/link_epi_to_expr")

gene<-strsplit(gsub("[a-z,/,.]","\\1",Sys.glob("temp/predictions/*"))[1],'_')[[1]][1]
i=1
list_positions<-c()
for (file in gsub("[a-z,/,.]","\\1",Sys.glob("temp/predictions/*"))){
    list_positions[i]<-strsplit(gsub("[a-z,/,.]","\\1",Sys.glob("temp/predictions/*")),'_')[[i]][2]
    i<-i+1
}

expression <- read.table(paste("temp/expression/",gene,".tsv",sep=''), head=T)

for (position in list_positions){
    predictions <- read.table(paste("temp/predictions/predictionsout",gene,"_",position,".fa.tsv",sep=''), head=T,row.names=1)
    #create correlation file
    if (position==list_positions[1])
        correlations <- data.frame(row.names=colnames(predictions)[5:dim(predictions)[2]])
    correlations<-cbind(correlations,data.frame(g=rep(0,919)))
    colnames(correlations)[colnames(correlations)=='g']<-paste("p",position,sep='')
    #add a column for expression
    predictions$expression<-predictions$metadata.ranges.id
    a=1
    for (i in row.names(predictions)){
        predictions$expression[a]<-expression[1,match(c(substring(i,4,10)),colnames(expression))]
        a<-a+1
    }

    for (i in colnames(predictions[5:(dim(predictions)[2]-1)])){
        correlations[i,paste("p",position,sep='')]<-cor(eval(parse(text = paste("predictions$",i,sep=""))),predictions$expression)
        #print(paste(i,cor(eval(parse(text = paste("predictions$",i,sep=""))),predictions$expression)))
        #print()
    }
}
warnings()
write.csv(correlations,file="correlations.csv")
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
