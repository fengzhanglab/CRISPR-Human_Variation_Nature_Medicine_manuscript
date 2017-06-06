rm(list = ls(.GlobalEnv), envir = .GlobalEnv)

library("reshape2")
library("ggplot2")

root_path <- "/broad/zhanglabdata/DAS_sandbox/_genoGen"
fig_path <- "Reference/ftp.broadinstitute.org/161006_query_vrer_pam_v0_2_1"
fout_path <- paste(root_path,"/",fig_path,"/sandbox",sep = "")
fin <- paste(root_path,"/",fig_path,"/_search.vrer.pq.csv",sep = "")
print(fin)

gff3_cols <- c('chrID','src','ftype','lb','ub','dot1','str','dot2','info')
af_cols <- c('af1e-5','af1e-4','af1e-3','af1e-2','af1e-1','af1','an')
het_cols <- c('het1e-5','het1e-4','het1e-3','het1e-2','het1e-1','het1','an')
n_cols <- c('ac','an','ptype')
all_cols <- c(gff3_cols,af_cols[1:6],het_cols[1:6],n_cols)
print(all_cols)

data <- read.csv(file=fin, header=FALSE, sep=",")
colnames(data) <- all_cols
head(data)

data_del <- data[data$ptype == 'del',]
data_add <- data[data$ptype == 'add',]

af_data_del <- data_del[af_cols]
af_data_add <- data_add[af_cols]
het_data_del <- data_del[het_cols]
het_data_add <- data_add[het_cols]

########################################
#figure1a
########################################
af_sum <- colSums(af_data_del)
het_sum <- colSums(het_data_del)

##########
df <- data.frame(
  group = c('del','no_del'),
  value = c(sum(af_sum[1:6])/af_sum[7],(af_sum[7]-sum(af_sum[1:6]))/af_sum[7])
  )

print(df)

bp <- ggplot(df, aes(x="", y=value, fill=group))+
geom_bar(width = 1, stat = "identity")

pie <- bp + coord_polar("y", start=0)+
scale_fill_brewer(palette = "Set1")+
theme_classic()
ggsave(filename=paste(fout_path,"/_fig1a-1.pdf",sep = ""))
write.csv(df,paste(fout_path,"/_fig1a-1.csv",sep = ""))

##########
df <- data.frame(
  group = rev(af_cols[1:6]),
  value = rev(af_sum[1:6]/af_sum[7])
  )

print(df)

bp <- ggplot(df, aes(x="", y=value, fill=rev(factor(group))))+
geom_bar(width = 1, stat = "identity")

pie <- bp + coord_polar("y", start=0)+
scale_fill_brewer(palette = "Set1")+
theme_classic()
ggsave(filename=paste(fout_path,"/_fig1a-2.pdf",sep = ""))
write.csv(df,paste(fout_path,"/_fig1a-2.csv",sep = ""))

########################################
#figure1b
########################################
af_sum <- colSums(af_data_del[1:6])
het_sum <- colSums(het_data_del[1:6])

hom_ratio <- c()
for (i in 1:length(af_sum)) {
  hom_ratio[i] <- (1-(het_sum[i]/af_sum[i]))
}

df <- data.frame(
  group = rev(af_cols[1:6]),
  value = rev(hom_ratio)
  )

print(df)

bp <- ggplot(df, aes(x=factor(group), y=value))+
geom_bar(stat = 'identity')+
scale_y_continuous(limits=c(0, 1), breaks=seq(0,1,0.2))+
theme_classic()
ggsave(filename=paste(fout_path,"/_fig1b.pdf",sep = ""))
write.csv(df,paste(fout_path,"/_fig1b.csv",sep = ""))

########################################
#figure1x
########################################
df <- data.frame(
  an = 1:length(data_del$ac),
  ac = 1:length(data_del$ac)
  )
df$an <- data_del$an
df$ac <- data_del$ac

dfc <- df[sample(nrow(df), 20000), ]
p1 <- ggplot(dfc, aes(x = an, y = ac))+ 
geom_point(size = 2)+
scale_y_continuous(trans='log2')+
scale_x_continuous(trans='log2')
ggsave(filename=paste(fout_path,"/_fig1x-an_scatter.eps",sep = ""))

ggplot(data=df[df$ac>0,], aes(ac/an))+ 
  geom_histogram(aes(y =..density..)) 
  # geom_density()
ggsave(filename=paste(fout_path,"/_fig1x-an_hist.pdf",sep = "")) 

########################################
#figure1x
########################################
df <- data.frame(
  an = 1:length(data_del$ac),
  ac = 1:length(data_del$ac)
  )
df$an <- data_del$ub-data_del$lb
df$ac <- data_del$ac

dfc <- df[sample(nrow(df), 20000), ]
p1 <- ggplot(dfc, aes(x = an, y = ac))+ 
geom_point(size = 2)+
scale_y_continuous(trans='log2')+
scale_x_continuous(trans='log2')
ggsave(filename=paste(fout_path,"/_fig1x-elen_scatter.eps",sep = ""))

ggplot(data=df[df$ac>0,], aes(ac/an))+ 
  geom_histogram(aes(y =..density..)) 
  # geom_density() 
ggsave(filename=paste(fout_path,"/_fig1x-elen_hist.pdf",sep = ""))


