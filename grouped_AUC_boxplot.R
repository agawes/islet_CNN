### CNN performance boxplots
library(ggplot2)
library(reshape)
library(grid)

setwd("~/Desktop/islets CNN/1000nets_results/")
perf=read.table("1000nets.AUC_ROC_PR.chr2.txt",h=T, row.names=1)
perf_melt=melt(perf, )
### add feature category
perf_melt$group=rep("Promoter",nrow(perf_melt))
perf_melt[grepl("H3K27ac|H3K4me1|H2K4me1|LMR", perf_melt$variable),]$group="Enhancer"
perf_melt[grepl("ATAC|DNase|H2A.Z", perf_melt$variable),]$group="Open Chromatin"
perf_melt[grepl("H3K36me3|H3K79me2|H2K4me2", perf_melt$variable),]$group="Active"
perf_melt[grepl("FOXA2|MAFB|NKX2.2|NKX6.1|PDX1|CTCF", perf_melt$variable),]$group="TF"
perf_melt[grepl("H3K27me3|H3K9me", perf_melt$variable),]$group="Repressed"

perf_melt_AUC=subset(perf_melt, grepl("_AUC",variable))
perf_melt_AUC$variable=gsub("_AUC","", perf_melt_AUC$variable)
perf_melt_AUC$variable=gsub("_","", perf_melt_AUC$variable)
perf_melt_AUC$variable=gsub("Acinar","aci_", perf_melt_AUC$variable)
perf_melt_AUC$variable=gsub("Alpha","a_", perf_melt_AUC$variable)
perf_melt_AUC$variable=gsub("Beta","b_", perf_melt_AUC$variable)
perf_melt_AUC$variable=gsub("Exocrine","e_", perf_melt_AUC$variable)

### add spaces in the start, so that the labels have the same length - 10 characters
perf_melt_AUC$variable[perf_melt_AUC$variable=="a_ATAC"]="    a_ATAC"
perf_melt_AUC$variable[perf_melt_AUC$variable=="a_H3K4me3"]=" a_H3K4me3"
perf_melt_AUC$variable[perf_melt_AUC$variable=="aci_ATAC"]="  aci_ATAC"
perf_melt_AUC$variable[perf_melt_AUC$variable=="ATAC"]="      ATAC"
perf_melt_AUC$variable[perf_melt_AUC$variable=="b_ATAC"]="  b_ATAC"
perf_melt_AUC$variable[perf_melt_AUC$variable=="b_H3K4me3"]=" b_H3K4me3"
perf_melt_AUC$variable[perf_melt_AUC$variable=="CTCF"]="        CTCF"
perf_melt_AUC$variable[perf_melt_AUC$variable=="DNase"]="      DNase"
perf_melt_AUC$variable[perf_melt_AUC$variable=="e_H3K4me3"]=" e_H3K4me3"
perf_melt_AUC$variable[perf_melt_AUC$variable=="FOXA2"]="     FOXA2"
perf_melt_AUC$variable[perf_melt_AUC$variable=="H2A.Z"]="     H2A.Z"
perf_melt_AUC$variable[perf_melt_AUC$variable=="H2K4me1"]="   H2K4me1"
perf_melt_AUC$variable[perf_melt_AUC$variable=="H2K4me2"]="   H2K4me2"
perf_melt_AUC$variable[perf_melt_AUC$variable=="H3K27ac"]="   H3K27ac"
perf_melt_AUC$variable[perf_melt_AUC$variable=="H3K27me3"]="  H3K27me3"
perf_melt_AUC$variable[perf_melt_AUC$variable=="H3K36me3"]="  H3K36me3"
perf_melt_AUC$variable[perf_melt_AUC$variable=="H3K4me1"]="   H3K4me1"
perf_melt_AUC$variable[perf_melt_AUC$variable=="H3K79me2"]="  H3K79me2"
perf_melt_AUC$variable[perf_melt_AUC$variable=="H3K9ac"]="    H3K9ac"
perf_melt_AUC$variable[perf_melt_AUC$variable=="H3K9me3"]="   H3K9me3"
perf_melt_AUC$variable[perf_melt_AUC$variable=="LMR"]="        LMR"
perf_melt_AUC$variable[perf_melt_AUC$variable=="MAFB"]="      MAFB"
perf_melt_AUC$variable[perf_melt_AUC$variable=="NKX2.2"]="    NKX2.2"
perf_melt_AUC$variable[perf_melt_AUC$variable=="NKX6.1"]="    NKX6.1"
perf_melt_AUC$variable[perf_melt_AUC$variable=="PDX1"]="      PDX1"
perf_melt_AUC$variable[perf_melt_AUC$variable=="UMR"]="       UMR"


library(ggplot2)
text_format <- element_text(face = "bold", color="black",size = 16)

p1 <- ggplot(subset(perf_melt_AUC, group=="Promoter"), aes(x=variable, y=value)) + 
  geom_boxplot(fill="red") + theme_bw() +ylab("AUROC") + xlab("") + ylim(c(0.65,1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12),legend.position = "none") +
  facet_grid(. ~ group) + theme(strip.background =element_rect(fill="red"), 
                                strip.text.x = element_text(size = 12, colour = "white", face = "bold"))
p2 <- ggplot(subset(perf_melt_AUC, group=="Enhancer"), aes(x=variable, y=value)) + 
  geom_boxplot(fill="orange") + theme_bw() +ylab("AUROC") + xlab("") + ylim(c(0.65,1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12),legend.position = "none") +
  facet_grid(. ~ group) + theme(strip.background =element_rect(fill="orange"), 
                                strip.text.x = element_text(size = 12, colour = "white", face = "bold"))
p3 <- ggplot(subset(perf_melt_AUC, group=="Repressed"), aes(x=variable, y=value)) + 
  geom_boxplot(fill="grey") + theme_bw() +ylab("AUROC") + xlab("") + ylim(c(0.65,1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12),legend.position = "none") +
  facet_grid(. ~ group) + theme(strip.background =element_rect(fill="grey"), 
                                strip.text.x = element_text(size = 12, colour = "white", face = "bold"))
p4 <- ggplot(subset(perf_melt_AUC, group=="Open Chromatin"), aes(x=variable, y=value)) + 
  geom_boxplot(fill="yellow3") + theme_bw() +ylab("AUROC") + xlab("") + ylim(c(0.65,1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12),legend.position = "none") +
  facet_grid(. ~ group) + theme(strip.background =element_rect(fill="yellow3"), 
                                strip.text.x = element_text(size = 12, colour = "white", face = "bold"))
p5 <- ggplot(subset(perf_melt_AUC, group=="TF"), aes(x=variable, y=value)) + 
  geom_boxplot(fill="skyblue") + theme_bw() +ylab("AUROC") + xlab("") + ylim(c(0.65,1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12),legend.position = "none") +
  facet_grid(. ~ group) + theme(strip.background =element_rect(fill="skyblue"), 
                                strip.text.x = element_text(size = 12, colour = "white", face = "bold"))
p6 <- ggplot(subset(perf_melt_AUC, group=="Active"), aes(x=variable, y=value)) + 
  geom_boxplot(fill="green3") + theme_bw() +ylab("AUROC") + xlab("") + ylim(c(0.65,1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12),legend.position = "none") +
  facet_grid(. ~ group) + theme(strip.background =element_rect(fill="green3"), 
                                strip.text.x = element_text(size = 12, colour = "white", face = "bold"))

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

pdf("AUC_by_group.pdf", width=10)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 16))) # 3 rows, 5 columns
print(p1, vp = vplayout(1, 1:6))  # the big plot covers rows 1:2 and cols 1:3
print(p2, vp = vplayout(1, 7:11))
print(p3, vp = vplayout(1, 12:16))
print(p4, vp = vplayout(2, 1:6))
print(p5, vp = vplayout(2, 7:12))
print(p6, vp = vplayout(2, 13:16))
dev.off()


##### same for PR but add random line based on 



perf_melt_PR=subset(perf_melt, grepl("_PR",variable))
perf_melt_PR$variable=gsub("_PR","", perf_melt_PR$variable)
perf_melt_PR$variable=gsub("_","", perf_melt_PR$variable)
perf_melt_PR$variable=gsub("Acinar","aci_", perf_melt_PR$variable)
perf_melt_PR$variable=gsub("Alpha","a_", perf_melt_PR$variable)
perf_melt_PR$variable=gsub("Beta","b_", perf_melt_PR$variable)
perf_melt_PR$variable=gsub("Exocrine","e_", perf_melt_PR$variable)

### add spaces in the start, so that the labels have the same length - 10 characters
perf_melt_PR$variable[perf_melt_PR$variable=="a_ATAC"]="    a_ATAC"
perf_melt_PR$variable[perf_melt_PR$variable=="a_H3K4me3"]=" a_H3K4me3"
perf_melt_PR$variable[perf_melt_PR$variable=="aci_ATAC"]="  aci_ATAC"
perf_melt_PR$variable[perf_melt_PR$variable=="ATAC"]="      ATAC"
perf_melt_PR$variable[perf_melt_PR$variable=="b_ATAC"]="  b_ATAC"
perf_melt_PR$variable[perf_melt_PR$variable=="b_H3K4me3"]=" b_H3K4me3"
perf_melt_PR$variable[perf_melt_PR$variable=="CTCF"]="        CTCF"
perf_melt_PR$variable[perf_melt_PR$variable=="DNase"]="      DNase"
perf_melt_PR$variable[perf_melt_PR$variable=="e_H3K4me3"]=" e_H3K4me3"
perf_melt_PR$variable[perf_melt_PR$variable=="FOXA2"]="     FOXA2"
perf_melt_PR$variable[perf_melt_PR$variable=="H2A.Z"]="     H2A.Z"
perf_melt_PR$variable[perf_melt_PR$variable=="H2K4me1"]="   H2K4me1"
perf_melt_PR$variable[perf_melt_PR$variable=="H2K4me2"]="   H2K4me2"
perf_melt_PR$variable[perf_melt_PR$variable=="H3K27ac"]="   H3K27ac"
perf_melt_PR$variable[perf_melt_PR$variable=="H3K27me3"]="  H3K27me3"
perf_melt_PR$variable[perf_melt_PR$variable=="H3K36me3"]="  H3K36me3"
perf_melt_PR$variable[perf_melt_PR$variable=="H3K4me1"]="   H3K4me1"
perf_melt_PR$variable[perf_melt_PR$variable=="H3K79me2"]="  H3K79me2"
perf_melt_PR$variable[perf_melt_PR$variable=="H3K9ac"]="    H3K9ac"
perf_melt_PR$variable[perf_melt_PR$variable=="H3K9me3"]="   H3K9me3"
perf_melt_PR$variable[perf_melt_PR$variable=="LMR"]="        LMR"
perf_melt_PR$variable[perf_melt_PR$variable=="MAFB"]="      MAFB"
perf_melt_PR$variable[perf_melt_PR$variable=="NKX2.2"]="    NKX2.2"
perf_melt_PR$variable[perf_melt_PR$variable=="NKX6.1"]="    NKX6.1"
perf_melt_PR$variable[perf_melt_PR$variable=="PDX1"]="      PDX1"
perf_melt_PR$variable[perf_melt_PR$variable=="UMR"]="       UMR"

#### read in the activation file & calculate class imbalance
a=read.table("chr1_2_last.act.txt",sep="\t",h=T,row.names=1)
table=t(sapply(colnames(a), function(x) table(a[,x])))
class_imb=apply(table, 1, function(x) x[2]/sum(x) )
names(class_imb)[2]="aci_ATAC"
names(class_imb)[3]="a_ATAC"
names(class_imb)[4]="a_H3K27me3"
names(class_imb)[5]="a_H3K4me3"
names(class_imb)[6]="b_ATAC"
names(class_imb)[7]="b_H3K27me3"
names(class_imb)[8]="b_H3K4me3"
names(class_imb)[11]="e_H3K27me3"
names(class_imb)[12]="e_H3K4me3"

text_format <- element_text(face = "bold", color="black",size = 16)

p1 <- ggplot(subset(perf_melt_PR, group=="Promoter"), aes(x=variable, y=value)) + 
  geom_boxplot(fill="red") + theme_bw() +ylab("AUPRC") + xlab("") + 
  ylim(c(0,1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12),legend.position = "none") +
  facet_grid(. ~ group) + theme(strip.background =element_rect(fill="red"), 
                                strip.text.x = element_text(size = 12, colour = "white", face = "bold")) 
df=data.frame(x=c(1:6), y=as.numeric(class_imb[c("UMR","H3K9ac","a_H3K4me3","b_H3K4me3", "e_H3K4me3","H3K4me3")]))
p1=p1 + geom_point(data=df, aes(x=x, y=y), size=2, shape=1)

p2 <- ggplot(subset(perf_melt_PR, group=="Enhancer"), aes(x=variable, y=value)) + 
  geom_boxplot(fill="orange") + theme_bw() +ylab("AUPRC") + xlab("") + 
  ylim(c(0,1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12),legend.position = "none") +
  facet_grid(. ~ group) + theme(strip.background =element_rect(fill="orange"), 
                                strip.text.x = element_text(size = 12, colour = "white", face = "bold"))
df=data.frame(x=c(1:4), y=as.numeric(class_imb[c("LMR","H2K4me1","H3K27ac","H3K4me1")]))
p2=p2 + geom_point(data=df, aes(x=x, y=y), size=2, shape=1)

p3 <- ggplot(subset(perf_melt_PR, group=="Repressed"), aes(x=variable, y=value)) + 
  geom_boxplot(fill="grey") + theme_bw() +ylab("AUPRC") + xlab("") + 
  ylim(c(0,1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12),legend.position = "none") +
  facet_grid(. ~ group) + theme(strip.background =element_rect(fill="grey"), 
                                strip.text.x = element_text(size = 12, colour = "white", face = "bold"))
df=data.frame(x=c(1:5), y=as.numeric(class_imb[c("H3K9me3","H3K27me3","a_H3K27me3","b_H3K27me3","e_H3K27me3")]))
p3=p3 + geom_point(data=df, aes(x=x, y=y), size=2, shape=1)
p3

p4 <- ggplot(subset(perf_melt_PR, group=="Open Chromatin"), aes(x=variable, y=value)) + 
  geom_boxplot(fill="yellow3") + theme_bw() +ylab("AUPRC") + xlab("") + 
  ylim(c(0,1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12),legend.position = "none") +
  facet_grid(. ~ group) + theme(strip.background =element_rect(fill="yellow3"), 
                                strip.text.x = element_text(size = 12, colour = "white", face = "bold"))
df=data.frame(x=c(1:6), y=as.numeric(class_imb[c("ATAC","DNase","H2A.Z","a_ATAC","aci_ATAC","b_ATAC")]))
p4=p4 + geom_point(data=df, aes(x=x, y=y), size=2, shape=1)
p4

p5 <- ggplot(subset(perf_melt_PR, group=="TF"), aes(x=variable, y=value)) + 
  geom_boxplot(fill="skyblue") + theme_bw() +ylab("AUPRC") + xlab("") +
  ylim(c(0,1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12),legend.position = "none") +
  facet_grid(. ~ group) + theme(strip.background =element_rect(fill="skyblue"), 
                                strip.text.x = element_text(size = 12, colour = "white", face = "bold"))
df=data.frame(x=c(1:6), y=as.numeric(class_imb[c("CTCF","MAFB","PDX1","FOXA2","NKX2.2","NKX6.1")]))
p5=p5 + geom_point(data=df, aes(x=x, y=y), size=2, shape=1)
p5

p6 <- ggplot(subset(perf_melt_PR, group=="Active"), aes(x=variable, y=value)) + 
  geom_boxplot(fill="green3") + theme_bw() +ylab("AUPRC") + xlab("") + 
  ylim(c(0,1)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12),legend.position = "none") +
  facet_grid(. ~ group) + theme(strip.background =element_rect(fill="green3"), 
                                strip.text.x = element_text(size = 12, colour = "white", face = "bold"))
df=data.frame(x=c(1:3), y=as.numeric(class_imb[c("H2K4me2","H3K36me3","H3K79me2")]))
p6=p6 + geom_point(data=df, aes(x=x, y=y), size=2, shape=1)
p6
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

pdf("PR_by_group.pdf", width=10)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 16))) # 3 rows, 5 columns
print(p1, vp = vplayout(1, 1:6))  # the big plot covers rows 1:2 and cols 1:3
print(p2, vp = vplayout(1, 7:11))
print(p3, vp = vplayout(1, 12:16))
print(p4, vp = vplayout(2, 1:6))
print(p5, vp = vplayout(2, 7:12))
print(p6, vp = vplayout(2, 13:16))
dev.off()

