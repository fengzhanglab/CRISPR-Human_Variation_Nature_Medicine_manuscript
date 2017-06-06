rm(list = ls(.GlobalEnv), envir = .GlobalEnv)

library("stringr")
library("reshape2")
library("ggplot2")
library("plyr")

fin <- "../_fig3a.csv"
print(fin)

data <- read.csv(file=fin, header=FALSE, sep=",")
data <- data[data$V8!="all",]
print(data)

p <- ggplot(data, aes(factor(V6), V5))+   
geom_bar(aes(fill = V8), position = "dodge", stat="identity")+
scale_fill_brewer(palette = "Set1")+
theme_classic()
# theme(axis.text.x=element_text(angle=90, hjust=1))
# ggsave(filename=paste(fout_path,"/_fig3a.pdf",sep = ""), width = 9, height = 3)
ggsave(filename="_fig3b_linear.pdf", width = 10, height = 5)
