rm(list = ls(.GlobalEnv), envir = .GlobalEnv)

library("stringr")
library("reshape2")
library("ggplot2")
library("plyr")

args<-commandArgs(TRUE)
fin <- c(args[1]);
fout <- c(args[2]);

data <- read.csv(file=fin, header=FALSE, sep=",")
print(data)

data <- data[data[,length(data[1,])]=='ALL',]

notplatinum <- data$X1>=1e-4
platinum <- data$X1<1e-4

df <- data.frame(ID = c('>=1e-4','<1e-4','<1e-4/total'))
gene <- 'CEP290'
temp <- data[grep(gene,data$V2),]
temp1 <- length(temp$V1[temp$V1>=1e-4])
temp2 <- length(temp$V1[temp$V1<1e-4])
temp3 <- temp2/(temp1+temp2)
df$CEP290 <- c(temp1,temp2,temp3)

gene <- 'CFTR'
temp <- data[grep(gene,data$V2),]
temp1 <- length(temp$V1[temp$V1>=1e-4])
temp2 <- length(temp$V1[temp$V1<1e-4])
temp3 <- temp2/(temp1+temp2)
df$CFTR <- c(temp1,temp2,temp3)

gene <- 'DMD'
temp <- data[grep(gene,data$V2),]
temp1 <- length(temp$V1[temp$V1>=1e-4])
temp2 <- length(temp$V1[temp$V1<1e-4])
temp3 <- temp2/(temp1+temp2)
df$DMD <- c(temp1,temp2,temp3)

gene <- 'G6PC'
temp <- data[grep(gene,data$V2),]
temp1 <- length(temp$V1[temp$V1>=1e-4])
temp2 <- length(temp$V1[temp$V1<1e-4])
temp3 <- temp2/(temp1+temp2)
df$G6PC <- c(temp1,temp2,temp3)

gene <- 'HBB'
temp <- data[grep(gene,data$V2),]
temp1 <- length(temp$V1[temp$V1>=1e-4])
temp2 <- length(temp$V1[temp$V1<1e-4])
temp3 <- temp2/(temp1+temp2)
df$HBB <- c(temp1,temp2,temp3)

gene <- 'IDUA'
temp <- data[grep(gene,data$V2),]
temp1 <- length(temp$V1[temp$V1>=1e-4])
temp2 <- length(temp$V1[temp$V1<1e-4])
temp3 <- temp2/(temp1+temp2)
df$IDUA <- c(temp1,temp2,temp3)

gene <- 'IL2RG'
temp <- data[grep(gene,data$V2),]
temp1 <- length(temp$V1[temp$V1>=1e-4])
temp2 <- length(temp$V1[temp$V1<1e-4])
temp3 <- temp2/(temp1+temp2)
df$IL2RG <- c(temp1,temp2,temp3)

gene <- 'PCSK9'
temp <- data[grep(gene,data$V2),]
temp1 <- length(temp$V1[temp$V1>=1e-4])
temp2 <- length(temp$V1[temp$V1<1e-4])
temp3 <- temp2/(temp1+temp2)
df$PCSK9 <- c(temp1,temp2,temp3)

gene <- 'PDCD1'
temp <- data[grep(gene,data$V2),]
temp1 <- length(temp$V1[temp$V1>=1e-4])
temp2 <- length(temp$V1[temp$V1<1e-4])
temp3 <- temp2/(temp1+temp2)
df$PDCD1 <- c(temp1,temp2,temp3)

gene <- 'SERPINA1'
temp <- data[grep(gene,data$V2),]
temp1 <- length(temp$V1[temp$V1>=1e-4])
temp2 <- length(temp$V1[temp$V1<1e-4])
temp3 <- temp2/(temp1+temp2)
df$SERPINA1 <- c(temp1,temp2,temp3)

gene <- 'TTR'
temp <- data[grep(gene,data$V2),]
temp1 <- length(temp$V1[temp$V1>=1e-4])
temp2 <- length(temp$V1[temp$V1<1e-4])
temp3 <- temp2/(temp1+temp2)
df$TTR <- c(temp1,temp2,temp3)

gene <- 'VEGFA'
temp <- data[grep(gene,data$V2),]
temp1 <- length(temp$V1[temp$V1>=1e-4])
temp2 <- length(temp$V1[temp$V1<1e-4])
temp3 <- temp2/(temp1+temp2)
df$VEGFA <- c(temp1,temp2,temp3)

temp1 <- length(data$V1[data$V1>=1e-4])
temp2 <- length(data$V1[data$V1<1e-4])
temp3 <- temp2/(temp1+temp2)
df$total <- c(temp1,temp2,temp3)

print(df)
write.csv(df,paste(gsub(".pdf",".csv",fout),sep = ""))

p <- ggplot(data,aes(x = factor(V2), y = V1))+ 
geom_hline(yintercept=c(1e-6,1e-5,1e-4,1e-3,1e-2,1e-1), colour="grey80", size = 0.2)+
geom_boxplot(position = "dodge")+
scale_colour_brewer(palette = "Set1")+
theme_classic()+
scale_y_log10(limits = c(1e-6, 1))
ggsave(filename=gsub(".pdf",".eps",fout), width = 9, height = 3)

data$V6 <- data$V1>0
data$V7 <- data$V1>=0
# print(data)

datac <- melt(rowsum(as.numeric(data$V6),data$V2)/rowsum(as.numeric(data$V7),data$V2))
# print(datac)

p <- ggplot(datac,aes(x = factor(Var1), y = value, stat = value))+ 
# p <- ggplot(data,aes(x = factor(V4), y = V7))+ 
geom_hline(yintercept=c(0.7,0.8,0.9,1), colour="grey80", size = 0.2)+
scale_y_continuous(limits=c(0.695, 1), breaks=seq(0.70,1,0.1))+
geom_point(size=5)+
theme_classic()
# scale_y_continuous(trans='log10')
ggsave(filename=gsub(".pdf","-3.eps",fout), width = 9, height = 1.5)


