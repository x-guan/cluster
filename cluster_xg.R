### This R code group RNA-Seq contigs into clusters based on their expression level
### Written by Xin Guan (github.com/x-guan)
library(tidyverse)
library(MBCluster.Seq)

rm(list=ls())
setwd("/Users/Guan/Xin/R/Github/cluster")

### generate demo data mc1:mc1000
dat_demo <- read.csv("/Users/Guan/Xin/Unix/SCG4_server/mca//180113_mca_getorf/mc_count_norm.csv", header=T)
dat_demo <- dat_demo[1:1000,]
write.csv(dat_demo, "mc_count_demo.csv", row.names=F)

### input RNA-Seq data and experimental design info
dat_1 <- read.csv("mc_count_demo.csv", header = T)
dat_2 <- dat_1[,4:27]   # RNA.table
contig_id <- dat_1[,1]   # GeneID
treatment <- rep(1:24, each=1)  # Treatment
# standardized RNA-Seq data
mydata <- RNASeq.Data(dat_2, Normalizer = NULL, treatment, contig_id)

####################################################################################################
### calculate AIC that determines the number of clusters
calcAIC <- function(i){
  set.seed(94305)
  nk <- i
  # monitor progess
  cat(paste("number of cluster = ", nk, "\n", sep = ""))
  c0 <- KmeansPlus.RNASeq(mydata,nK=nk)$centers 
  # use EM algorithm to cluster genes; cls gives the group id of each gene
  cls <- Cluster.RNASeq(data = mydata, model = "nbinom", centers = c0, method = "EM")
  # calculate log(likelihood)
  lglk <- lglk.cluster(mydata,model="nbinom",cluster=cls$cluster)
  nG <- length(contig_id)
  I <- length(unique(treatment))
  np <- nG*(nk + 1) + nk*I -1
  # calculate AIC
  AIC <- -2*(lglk*nG - np)
  cat(paste("AIC = ", AIC, "\n", sep = ""))
  return(list(AIC=AIC,cls=cls))
}
# execute calcAIC
nk_vec <- c(2:100)
AIC_re = vector("list", length(nk_vec))
for(i in 1:length(nk_vec)){
  AIC_re[[i]] <- calcAIC(nk_vec[i])
  }
AIC_re_2 <- data.frame(nK = nk_vec, AIC = sapply(AIC_re, function(x)x$AIC))
write.csv(AIC_re_2, "cl_1_aic.csv", row.names=F)
AIC_re_2 <- read.csv("cl_1_aic.csv", header=T)
ggplot(AIC_re_2, aes(x=nK, y=AIC)) +
  theme_bw() +
  geom_point(position = position_dodge(width = 0.2)) +
  # stat_smooth(method=lm, formula=y~poly(x,2), se=F, size=0.2) +
  ggtitle("AIC plot") +
  # ylab("AIC") +
  # xlab('Number of clusters') +
  scale_y_continuous("AIC") +
  scale_x_continuous("Number of clusters", breaks=seq(10,100, by=10)) +
  ggsave("cl_1_aic_plot.pdf", height=4, width=6)
# the smallest AIC value (380131) come out at 47 clusters

####################################################################################################
### asign clusters
set.seed(94305)
# choose the number of cluster centers to initialize the clustering
nk <- 47
c0 <- KmeansPlus.RNASeq(mydata,nK=nk)$centers 
# use EM algorithm to cluster genes
cls <- Cluster.RNASeq(data=mydata, model="nbinom", centers=c0, method="EM")$cluster
# give cls info to each contig
dat_3 <- cbind(dat_1, cls)
write.csv(dat_3, "cl_2_info.csv", row.names=F)
# bulild a tree structure
tree_info <- Hybrid.Tree(data=mydata, cluste=cls, model="nbinom")
write.csv(tree_info, file = "cl_3_tree.csv", row.names=F)
# plot the tree structure
pdf("cl_3_tree.pdf", height=9, width=12)
plotHybrid.Tree(merge=tree_info, cluster=cls, logFC=mydata$logFC, tree.title=NULL, colorful=T)
dev.off()

####################################################################################################
### plot contigs in each cluster
cls_info <- read.csv("cl_2_info.csv", header=T)
# check cluster info
# cls_info[cls_info$contig_id == "mc10881",]$cls
cls_info[,4:27] <- t(apply(cls_info[,4:27]+1, 1, function(x){log2(x/mean(x))}))
# plot
cls_info_2 <- cls_info %>%
  select(-c(contig, length)) %>%
  gather(organ, relative_exp, -c(contig_id, cls)) %>%
  mutate(organ = factor(organ, levels = c(paste0("l",1:12), paste0("p",1:3), paste0("s",1:3),
                                          paste0("rh",1:3), paste0("r",1:3)))) 
ymax <- max(cls_info_2$relative_exp)
ymin <- min(cls_info_2$relative_exp)
ggplot(cls_info_2, aes(x=organ, y=relative_exp, group=contig_id)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6)) +
  theme(axis.text.y=element_text(size=7)) +
  geom_rect(aes(xmin="l1", xmax="l6", ymin=-Inf, ymax=Inf), fill="grey", alpha=0.1) +
  geom_rect(aes(xmin="l7", xmax="l12", ymin=-Inf, ymax=Inf), fill="grey", alpha=0.1) +
  geom_rect(aes(xmin="p1", xmax="p3", ymin=-Inf, ymax=Inf), fill="grey", alpha=0.1) +
  geom_rect(aes(xmin="s1", xmax="s3", ymin=-Inf, ymax=Inf), fill="grey", alpha=0.1) +
  geom_rect(aes(xmin="rh1", xmax="rh3", ymin=-Inf, ymax=Inf), fill="grey", alpha=0.1) +
  geom_rect(aes(xmin="r1", xmax="r3", ymin=-Inf, ymax=Inf), fill="grey", alpha=0.1) +
  geom_line(color="black", size=0.1) +
  stat_summary(fun.y=mean, color="darkred", size=1, geom="line", group=1) +
  ggtitle("Clusters of gene expression") +
  scale_y_continuous("Relative expression level", limits=c(ymin, ymax), breaks=seq(-16, 8, by=2)) +
  xlab('Menispermum organs') +
  facet_wrap(~cls, ncol=10, drop=T) +
  ggsave("cl_4_cluster.pdf", height=14, width=20)







