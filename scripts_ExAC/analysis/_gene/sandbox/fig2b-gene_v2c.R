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

data$V6 <- data$V1+data$V2+data$V3

#filtering
data <- data[data$V5!=100,]
data <- data[data$V5<1e5,]

print(data)

# data.m = melt(data,id=c("V2","V3"))
# print(data.m)

# p <- ggplot(data,aes(x = V3, y = V1, group = factor(V2), color = factor(V2)))+ 
# geom_line(size = 1)+
# scale_colour_brewer(palette = "Paired")+
# # scale_y_continuous(limits=c(0, 0.15), breaks=seq(0,0.15,0.02))+
# theme_classic()
# # scale_y_continuous(trans='log10')
# # dev.new(width=9, height=3)
# ggsave(filename=fout, width = 7, height = 7)

p <- ggplot(data,aes(x = V5, y = V6, group = factor(V4), color = factor(V4)))+ 
geom_line(size = 1)+
scale_colour_brewer(palette = "Paired")+
# scale_y_continuous(limits=c(0, 0.15), breaks=seq(0,0.15,0.02))+
theme_classic()
# scale_x_continuous(trans='log10')+
# scale_y_continuous(trans='log10')
# dev.new(width=9, height=3)
ggsave(filename=fout, width = 7, height = 7)

genes <- unique(data$V4)
for (i in 1:length(genes)) {
	data.m <- melt(data[data$V4==genes[i],],id=c("V4","V5"))
	data.m$variable <- factor(data.m$variable, levels=c("V3","V2","V1","V6"))
	p <- ggplot(data.m,aes(x = V5, y = value, group = factor(variable), color = factor(variable)))+ 
	geom_line(size = 1)+
	scale_colour_brewer(palette = "Set1")+
	scale_y_continuous(limits=c(0, 0.645), breaks=seq(0,0.645,0.1))+
	theme_classic()
	ggsave(filename=gsub('.pdf',paste('_',genes[i],'.pdf',sep=''),fout), width = 7, height = 7)
}


