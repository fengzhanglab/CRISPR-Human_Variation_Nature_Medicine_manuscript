rm(list = ls(.GlobalEnv), envir = .GlobalEnv)

library("stringr")
library("reshape2")
library("ggplot2")
library("plyr")

args<-commandArgs(TRUE)
fin <- c(args[1]);
print(fin)


data <- read.csv(file=fin, header=FALSE, sep=",")
head(data)

data.m = melt(data,id=c("V27","V28","V29"))

data.m <- ddply(data.m, .(variable), transform, rescale=value)
print(data.m)

p <- ggplot(data.m,aes(x = variable, y = value, group = factor(V29), color = rev(factor(V29))))+ 
geom_line(size = 1)+
scale_y_log10(limits = c(3e-4, 0.1))+
# scale_y_continuous(trans='log10')+
scale_colour_brewer(palette = "Set1")+
theme_classic()
# dev.new(width=9, height=3)
ggsave(filename="_fig2a-v2_log10.pdf", width = 9, height = 3)

p <- ggplot(data.m,aes(x = variable, y = value, group = factor(V29), color = rev(factor(V29))))+ 
geom_line(size = 1)+
scale_colour_brewer(palette = "Set1")+
theme_classic()
# dev.new(width=9, height=3)
ggsave(filename="_fig2a-v2_linear.pdf", width = 9, height = 3)


