rm(list = ls(.GlobalEnv), envir = .GlobalEnv)

library("stringr")
library("reshape2")
library("ggplot2")
library("plyr")

args<-commandArgs(TRUE)
fin <- c(args[1]);
ain <- c(args[2]);
fout <- c(args[3]);

data <- read.csv(file=fin, header=FALSE, sep=",")
data$V7 <- as.numeric(as.character(data$V7))

annot <- read.csv(file=ain, header=FALSE, sep=",")
# annot$V4 <- as.numeric(as.character(annot$V4))
# annot$V5 <- as.numeric(as.character(annot$V5))

for (i in 1:length(annot$V1)) {
	print(i)
	lb <- annot$V4[i]
	ub <- annot$V5[i]
	elen <- ub-lb

	datam <- data[data$V7>=lb,]
	datam <- datam[datam$V7<=ub,]
	if (length(datam$V1) == 0) {
		print('no values for annot row:')
		print(annot[i,])
	} else {

		if (annot$V6[i] == '+') {
			datam <- datam[order(datam$V7),]
			lb2 <- lb
			ub2 <- ub
		} else if (annot$V6[i] == '-') {
			datam <- datam[order(-datam$V7),]
			lb2 <- ub
			ub2 <- lb
		}

		datam <- melt(datam,id=c("V5","V6","V7","V8","V9","V10","V11"))	

		p <- ggplot(datam,aes(x = V7, y = value, group = factor(variable), colour = factor(variable)))+ 
		geom_hline(yintercept=c(1,10,100), colour="grey80", size=0.2)+
		geom_point(size=3, alpha=0.8)+
		scale_colour_brewer(palette="Set1")+
		scale_y_log10(limits=c(1, 1000))+
		theme_classic()+
		theme(axis.text.x=element_text(angle=90, hjust=1))+
		geom_vline(xintercept=lb:ub, colour="grey80", linetype="longdash", size=0.1)

		width <- 12*(elen/100)
		if (width > 48) {width <- 48}
		print(gsub(".pdf", paste("_",toString(i),".pdf",sep=''), fout))
		ggsave(filename=gsub(".pdf", paste("_",toString(i),"-mm_v3.pdf",sep=''), fout), width = width , height = 3)	
	}	
}








