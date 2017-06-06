rm(list = ls(.GlobalEnv), envir = .GlobalEnv)

library("stringr")
library("reshape2")
library("ggplot2")
library("plyr")

args<-commandArgs(TRUE)
fin <- c(args[1]);
fout <- c(args[2]);

data <- read.csv(file=fin, header=FALSE, sep=",")
head(data)

data.m = melt(data,id=c("V24","V25"))
print(data.m)

p <- ggplot(data.m,aes(x = variable, y = value, group = factor(V25), color = rev(factor(V25))))+ 
geom_line(size = 1)+
scale_colour_brewer(palette = "Set1")+
scale_y_continuous(limits=c(0, 0.30), breaks=seq(0,0.30,0.02))+
theme_classic()
# scale_y_continuous(trans='log10')
# dev.new(width=9, height=3)
ggsave(filename=fout, width = 9, height = 3)

# p <- ggplot(data.m,aes(x = variable, y = value, group = factor(V25), color = rev(factor(V25))))+ 
# geom_line(size = 1)+
# scale_colour_brewer(palette = "Set1")+
# # scale_y_continuous(limits=c(0, 0.12), breaks=seq(0,0.12,0.02))+
# theme_classic()+
# scale_y_continuous(trans='log10')
# # dev.new(width=9, height=3)
# ggsave(filename=fout, width = 9, height = 3)


