rm(list = ls(.GlobalEnv), envir = .GlobalEnv)

library("stringr")
library("reshape2")
library("ggplot2")
library("plyr")

root_path <- "/broad/zhanglabdata/DAS_sandbox/_genoGen"
fig_path <- "Reference/ftp.broadinstitute.org/160924_query_spcas9NGG_gene_v0_2_1"
fin <- paste(root_path,"/",fig_path,"/_search-PCSK9_ENSE00001279167.spcas9.gq.csv",sep = "")
fout_path <- paste(root_path,"/",fig_path,"/sandbox",sep = "")
print(fin)

data <- read.csv(file=fin, header=FALSE, sep=",")
print(data)

print(length(data[1,]))
print(data[,length(data[1,])-2]=='ALL')

data <- data[data[,length(data[1,])-2]=='ALL',]
data$V41 <- paste(data$V33,data$V34,sep="-")
head(data)
tail(data)

# for (i in 10:32) {
# 	data[data[,i]<0.000001,i] <- 7
# 	data[data[,i]<0.00001,i] <- 6
# 	data[data[,i]<0.0001,i] <- 5
# 	data[data[,i]<0.001,i] <- 4
# 	data[data[,i]<0.01,i] <- 3
# 	data[data[,i]<0.1,i] <- 2
# 	data[data[,i]<=1,i] <- 1 
# }

data <- data[,c(10:32,41)]
head(data)
tail(data)

data.m = melt(data, id=c("V41"))
data.m <- rev(data.m[order(data.m$V41),])
print(data.m)

# p <- ggplot(data.m[data.m$value > 0,],aes(variable,factor(V33)))+ 
p <- ggplot(data.m,aes(variable,factor(V41)))+ 
geom_tile(aes(fill=value), colour = "white")+ 
theme_classic()+
scale_fill_gradient(low="dodgerblue2",high="firebrick2",na.value="grey87",trans='log10')+
scale_y_discrete(limits = rev(levels(factor(data.m$V41))))
# dev.new(width=9, height=3)
ggsave(filename=paste(fout_path,"/_fig2a_log10.pdf",sep = ""), width = 6, height = 10)


