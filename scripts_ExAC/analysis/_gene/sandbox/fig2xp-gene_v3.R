rm(list = ls(.GlobalEnv), envir = .GlobalEnv)

library("stringr")
library("reshape2")
library("ggplot2")
library("plyr")

args<-commandArgs(TRUE)
fin <- c(args[1]);
fout <- c(args[2]);

data <- read.csv(file=fin, header=FALSE, sep=",")

data <- data[data[,length(data[1,])]=='ALL',]
# data$V1[data$V1<1e-6] <- 1e-6
# data$V6 <- data$V1 >= 1e-4

print(data)

# data$V3 <- as.numeric(as.character(data$V3))
# data <- data[order(data$V3),]

data$V5 <- paste(data$V3,data$V4,sep="-")
data$V5 <- factor(data$V5, levels=unique(data$V5))

bin_size <- 30
n_bins <- as.integer(ceiling(length(data$V5)/bin_size))
print (n_bins)
for (i in 1:n_bins) {
	print(i)
	lb <- ((bin_size*(i-1))+1)
	ub <- (bin_size*i)
	if (lb < 1) {lb <- 1}
	if (ub > length(data$V5)) {ub <- length(data$V5)}

	current <- data[lb:ub,]
	# current$V6 <- factor(current$V6, levels=c("TRUE","FALSE"))

	p <- ggplot(current,aes(x = V5, y = V1))+ 
	geom_hline(yintercept=c(1e-6,1e-5,1e-3,1e-2,1e-1), colour="grey80", size = 0.2)+
	geom_hline(yintercept=c(1e-4), colour="grey60", size = 0.5)+
	geom_point(size = 3, alpha = 0.8)+
	scale_colour_brewer(palette = "Set1")+
	scale_y_log10(limits = c(1e-6, 1))+
	theme_classic()+
	theme(axis.text.x=element_text(angle=60, hjust=1))
	# geom_vline(xintercept = lb:ub, colour="grey50", linetype = "dash", size = 0.5)

	print(gsub(".pdf", paste("_",toString(i),".pdf",sep=''), fout))
	ggsave(filename=gsub(".pdf", paste("_",toString(i),".pdf",sep=''), fout), width = 8, height = 3)	
}


