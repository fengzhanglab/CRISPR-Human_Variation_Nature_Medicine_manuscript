rm(list = ls(.GlobalEnv), envir = .GlobalEnv)

library("stringr")
library("reshape2")
library("ggplot2")
library("plyr")

fin <- "../_fig3c.csv"
print(fin)

data <- read.csv(file=fin, header=FALSE, sep=",")
print(data)

p <- ggplot(data, aes(factor(V6), V5))+   
geom_bar(aes(fill = rev(factor(V8))), position = "dodge", stat="identity")+
theme_classic()+
scale_fill_brewer(palette = "Set1")
# theme(axis.text.x=element_text(angle=90, hjust=1))
# ggsave(filename=paste(fout_path,"/_fig3a.pdf",sep = ""), width = 9, height = 3)
ggsave(filename="_fig3c_linear.pdf", width = 10, height = 5)
