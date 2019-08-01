## first get the number of annotated motifs in each net:
## for I in fs*_runs/TF_motifs/iter*/table.txt; do printf "$I\t"; awk '{if(NR>1 && $3 != "."){print $0}}' $I | wc -l; done > annotated_motifs.count_per_net.txt

### count how many motifs in each net have std>0 - these are the motifs that help with predictions
## for I in fs*_runs/TF_motifs/iter*/table.txt; do printf "$I\t"; awk '{if(NR>1 && $6>0){print $0}}' $I | wc -l; done > motifs_with_std.count_per_net.txt

## count how many novel/unannotated motifs we discover
### for I in fs*_runs/TF_motifs/iter*/table.txt; do printf "$I\t"; awk '{if(NR>1 && $6>0 && $3 != "."){print $0}}' $I | wc -l; done > motifs_with_std.unannotated.count_per_net.txt

library(ggplot2)

setwd("~/Desktop/islets CNN/TF_summary/")
motifs=read.table("annotated_motifs.count_per_net.txt")
motifs$fs=gsub("fs","",gsub("_runs.+","",motifs$V1))
motifs$fs=as.numeric(motifs$fs)

### look at the distributions
summary(motifs$V2)
aggregate(motifs$V2, by=list(fs=motifs$fs), FUN=mean)

## make boxplot of known motifs detected by fs
motifs$fs=factor(motifs$fs)

text_format <- element_text( color="black",size = 16)

p1=ggplot(motifs, aes(x=fs, y=V2)) + geom_boxplot(fill="grey") + theme_bw() + xlab("Filter size") + ylab("Annotated filters") +
  theme(axis.text.x=text_format,axis.text.y=text_format, axis.title=text_format) 

motifs_std=read.table("motifs_with_std.count_per_net.txt")
motifs_std$fs=gsub("fs","",gsub("_runs.+","",motifs_std$V1))
motifs_std$fs=as.numeric(motifs_std$fs)

### look at the distributions
summary(motifs_std$V2)
aggregate(motifs_std$V2, by=list(fs=motifs_std$fs), FUN=mean)

## make boxplot of known motifs detected by fs
motifs_std$fs=factor(motifs_std$fs)
p2=ggplot(motifs_std, aes(x=fs, y=V2)) + geom_boxplot(fill="grey") + theme_bw() + xlab("Filter size") + ylab("Filters influencing predictions (std>0)") +
  theme(axis.text.x=text_format,axis.text.y=text_format, axis.title=text_format) 

### check how many of these filters were unannotated to TF's
motifs_std_un=read.table("motifs_with_std.unannotated.count_per_net.txt")
motifs_std_un$fs=gsub("fs","",gsub("_runs.+","",motifs_std_un$V1))
motifs_std_un$fs=as.numeric(motifs_std_un$fs)

### look at the distributions
summary(motifs_std_un$V2)
aggregate(motifs_std_un$V2, by=list(fs=motifs_std_un$fs), FUN=mean)

## make boxplot of known motifs detected by fs
motifs_std_un$fs=factor(motifs_std_un$fs)
ggplot(motifs_std_un, aes(x=fs, y=V2)) + geom_boxplot() + theme_bw() + xlab("Filter size") + ylab("Filters influencing predictions (std>0)") +
  theme(axis.text.x=text_format,axis.text.y=text_format, axis.title=text_format) 

p3 = ggplot(motifs_std, aes(x=fs, y=V2)) + geom_boxplot(fill="grey") + theme_bw() + xlab("Filter size") + ylab("Filters influencing predictions (std>0)") +
  theme(axis.text.x=text_format,axis.text.y=text_format, axis.title=text_format) + ylim(0,250 )
p3 + geom_boxplot(data=motifs_std_un, aes(x=fs, y=V2)) + geom_boxplot(data=motifs, aes(x=fs, y=V2), fill="indianred") 


motifs$label="Annotated"
motifs_std$label="Informative"
motifs_std_un$label="Informative\nUnannotated"

all_motifs=rbind(motifs, motifs_std)
all_motifs=rbind(all_motifs, motifs_std_un)

pdf("Filters.summary.pdf", width=6)
p3 = ggplot(motifs_std, aes(x=fs, y=V2)) + geom_boxplot(fill="grey") + theme_bw() + xlab("Filter size") + ylab("CNN Filters") +
  theme(axis.text.x=text_format,axis.text.y=text_format, axis.title=text_format) + ylim(0,250 )
p3 + geom_boxplot(data=motifs_std_un, aes(x=fs, y=V2)) + geom_boxplot(data=motifs, aes(x=fs, y=V2), fill="indianred") 

ggplot(all_motifs, aes(x=fs, y=V2, by=label)) + geom_boxplot(aes(fill=label)) + theme_bw() + xlab("Filter size") + ylab("CNN Filters") +
  theme(axis.text.x=text_format,axis.text.y=text_format, axis.title=text_format, legend.text = text_format, legend.title = text_format) + 
  scale_fill_manual(values=c("indianred","grey","white"), breaks=c("Annotated","Informative","Informative\nUnannotated" ), name="Filters")
dev.off()