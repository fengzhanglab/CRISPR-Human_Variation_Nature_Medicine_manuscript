rm(list = ls(.GlobalEnv), envir = .GlobalEnv)

library("stringr")
library("reshape2")
library("ggplot2")
library("plyr")

args<-commandArgs(TRUE)
fin <- c(args[1]);
fout <- c(args[2]);

data <- read.csv(file=fin, header=FALSE, sep=",")
# data <- data[data$V4<3000,]
data$V7 <- as.numeric(as.character(data$V7))
print(cor(data$V1,data$V2))
print(cor(data$V2,data$V3))
print(cor(data$V3,data$V4))

p <- ggplot(data,aes(x = V3, y = V4, group(factor(V10)), colour(factor(V10))))+ 
geom_point(size=1, alpha=0.8)+
scale_colour_brewer(palette="Paired")+
# scale_x_log10(limits=c(1, 1000))+
# scale_y_log10(limits=c(1, 1000))+
scale_x_log10()+
scale_y_log10()+
theme_classic()

ggsave(filename=fout, width = 7 , height = 7)	








