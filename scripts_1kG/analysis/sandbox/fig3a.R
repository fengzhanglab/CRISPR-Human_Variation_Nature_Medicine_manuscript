rm(list = ls(.GlobalEnv), envir = .GlobalEnv)

library("stringr")
library("reshape2")
library("ggplot2")
library("plyr")

fin <- paste("../_fig3a.csv",sep = "")
print(fin)

data <- read.csv(file=fin, header=FALSE, sep=",")
data <- data[data$V8=="all",]
print(data)

data.m = melt(data,id=c("V5","V6","V7","V8"))
print(data.m)

p <- ggplot(data.m,aes(x = factor(V6), y = value, group = factor(variable), color = factor(variable)))+ 
geom_point(size=4,alpha=1)+
scale_y_continuous(trans='log10')+
# scale_color_hue(l=65, c=150)+
theme_classic()+
theme(axis.text.x=element_text(angle=90, hjust=1))+
scale_colour_brewer(palette = "Set1")

# ggsave(filename=paste(fout_path,"/_fig3a.pdf",sep = ""), width = 9, height = 3)
ggsave(filename="_fig3a_log.eps", width = 10, height = 4)

# p <- ggplot(data.m,aes(x = factor(V8), y = value, group = factor(variable), color = factor(variable)))+ 
# geom_point(size=12,alpha=0.8)+
# # scale_color_hue(l=65, c=150)+
# theme_classic()+
# theme(axis.text.x=element_text(angle=90, hjust=1))+
# scale_colour_brewer(palette = "Set1")
# # ggsave(filename=paste(fout_path,"/_fig3a.pdf",sep = ""), width = 9, height = 3)
# ggsave(filename=paste(fout_path,"/_fig3a_linear.pdf",sep = ""), width = 8, height = 5)
