rm(list = ls(.GlobalEnv), envir = .GlobalEnv)

library("stringr")
library("reshape2")
library("ggplot2")
library("plyr")

args<-commandArgs(TRUE)
fin <- c(args[1]);
fout <- c(args[2]);

data <- read.csv(file=fin, header=FALSE, sep=",")

data$V7 <- as.numeric(as.character(data$V7))
data <- data[order(data$V7),]
data$V12 <- paste(data$V7,data$V9,sep="-")
data$V12 <- factor(data$V12, levels=unique(data$V12))
print(data)

bin_size <- 90
n_bins <- as.integer(ceiling(length(data$V12)/bin_size))
print(n_bins)
for (i in 1:n_bins) {
	print(i)
	lb <- ((bin_size*(i-1))+1)
	ub <- (bin_size*i)
	if (ub > length(data$V12)) {ub <- length(data$V12)}

	print(data[lb:ub,])
	p <- ggplot(data[lb:ub,],aes(x = V12, y = V5))+ 
	geom_hline(yintercept=c(1,10,100), colour="grey80", size = 0.2)+
	geom_point(size = 3, alpha = 0.7)+
	# scale_colour_brewer(palette = "Set1")+
	scale_y_log10(limits = c(1, 1000))+
	theme_classic()+
	theme(axis.text.x=element_text(angle=90, hjust=1))
	geom_vline(xintercept = lb:ub, colour="grey50", linetype = "dash", size = 0.5)

	print(gsub(".pdf", paste("_",toString(i),".pdf",sep=''), fout))
	ggsave(filename=gsub(".pdf", paste("_",toString(i),".pdf",sep=''), fout), width = 12 , height = 3)	

	datam <- data[lb:ub,]
	datam <- melt(datam,id=c("V5","V6","V7","V8","V9","V10","V11","V12"))
	p <- ggplot(datam,aes(x = V12, y = value, group = factor(variable), colour = factor(variable)))+ 
	geom_hline(yintercept=c(1,10,100), colour="grey80", size = 0.2)+
	geom_point(size = 3, alpha = 0.7)+
	scale_colour_brewer(palette = "Set1")+
	scale_y_log10(limits = c(1, 1000))+
	theme_classic()+
	theme(axis.text.x=element_text(angle=90, hjust=1))
	geom_vline(xintercept = lb:ub, colour="grey50", linetype = "dash", size = 0.5)

	print(gsub(".pdf", paste("_",toString(i),".pdf",sep=''), fout))
	ggsave(filename=gsub(".pdf", paste("_",toString(i),"-mm.pdf",sep=''), fout), width = 12 , height = 3)		
}

# bin_size <- 200
# root <- data$V7[1]
# data$V13 <- 0
# for (i in 1:length(data$V12)) {

# }








