rm(list = ls(.GlobalEnv), envir = .GlobalEnv)

library("reshape2")
library("ggplot2")

args<-commandArgs(TRUE)
fin <- c(args[1]);
print(fin)

data <- read.csv(file=fin, header=FALSE, sep=",")
colnames(data) <- c('af_sum','het_sum','var_total','targ_total','bin')
head(data)

########################################
#figure1b
########################################
data$hom_ratio <- 1-(data$het_sum/data$af_sum)

print(data)

df <- data.frame(
  group = rev(data$bin),
  value = rev(data$hom_ratio)
  )

print(df)

bp <- ggplot(df, aes(x=factor(group), y=value))+
geom_bar(stat = 'identity')+
scale_y_continuous(limits=c(0, 1), breaks=seq(0,1,0.2))+
theme_classic()
ggsave(filename="_fig2a-2.pdf")


