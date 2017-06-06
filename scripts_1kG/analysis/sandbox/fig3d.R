rm(list = ls(.GlobalEnv), envir = .GlobalEnv)

library("stringr")
library("reshape2")
library("ggplot2")
library("plyr")

fin <- "../_fig3d.csv"
print(fin)

data <- read.csv(file=fin, header=FALSE, sep=",")
data <- data[data$V2507==4,]
# head(data)

data_hist <- data.frame()
for (i in 1:length(data[,1])) {
	lb <- round((min(data[i,1:2504])-.1)*10)/10
	ub <- round((max(data[i,1:2504])+.1)*10)/10
	hist_temp <- hist(data.matrix(data[i,1:2504]),breaks=seq(lb,ub,by=0.1))
	# print(hist_temp)

	print(length(hist_temp$counts))
	print(length(hist_temp$density))
	print(length(hist_temp$mids))

	data_temp <- data.frame(
	  counts = hist_temp$counts,
	  density = hist_temp$density,
	  mids = hist_temp$mids
	  )
	data_temp$id <- data[i,2505]

	data_hist <- rbind(data_hist,data_temp)
}
print(data_hist)

# data.m <- melt(data)
# head(data.m)

p <- ggplot(data_hist, aes(factor(id), mids, color = counts))+   
geom_point(size = 10)+
scale_colour_gradient(high="firebrick1", low="dodgerblue1")+
theme_classic()
# theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave(filename="_fig3d_linear_v1.eps", width = 10, height = 5)

p <- ggplot(data_hist, aes(factor(id), mids, color = counts))+   
geom_point(size = 10)+
scale_colour_gradient(high="firebrick2", low="dodgerblue2")+
theme_classic()
# theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave(filename="_fig3d_linear_v2.eps", width = 10, height = 5)

p <- ggplot(data_hist, aes(factor(id), mids, size = counts))+   
geom_point()+
theme_classic()
# theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave(filename="_fig3d_linear_v3.eps", width = 10, height = 6)


