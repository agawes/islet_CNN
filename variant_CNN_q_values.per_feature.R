library(data.table)
library(pbapply)
library(Biobase)

res_dir="/well/got2d/agata/Basset/Full_islet_dataset/"
files_0.1<-list.files(path=res_dir, pattern=".HRC_below_0.1.snp_predict.txt.ann.pred_diff$",recursive=TRUE)
files<-list.files(path=res_dir, pattern=".HRC_0.1_0.5.snp_predict.txt.ann.pred_diff$",recursive=TRUE)
files_0.5<-list.files(path=res_dir, pattern=".HRC_above_0.5.snp_predict.txt.ann.pred_diff$",recursive=TRUE)

features=read.table("data/samples.txt")

### read in all data
cnn_res = list()
for(f in 1:1000){
    name = gsub("/HRC_credset_predictions/iter","",gsub(".HRC_below_0.1.snp_predict.txt.ann.pred_diff","",files_0.1[f]))
    df = fread(paste0(res_dir,files_0.1[f]),sep="\t")
    df2=fread(paste0(res_dir,files[f]),sep="\t")
    df=rbind(df,df2)
    df3=fread(paste0(res_dir,files_0.5[f]),sep="\t")
	df=rbind(df,df3)
	df=df[!duplicated(df$V1),]
	names(df) = c("SNP_allele",as.character(features$V1))
	cnn_res[[name]] = data.table(df)
	print(f)
}

### for each variant:
# 	- calculate the score-diff 
#	- calculate mean score for each feature across all networks
#	- calculate rank product

df_all_features = data.frame(variant=cnn_res[[1]]$SNP_allele)

for (feature in colnames(cnn_res[[1]])[2:31] ){
	feature_scores = subListExtract(cnn_res, feature)
	feature_scores_m = do.call(rbind, feature_scores)
	feature_mean = pbapply(feature_scores_m, 2, mean)

	feature_rank = pbapply(feature_scores_m, 1, function(x) frank(1-x))
	geo_mean_rank = pbapply(feature_rank, 1, function(x) exp(mean(log(x))))
	
	df_all_features = cbind(df_all_features, feature_mean, geo_mean_rank)
	names(df_all_features)[ncol(df_all_features)-1]=paste0(feature,"_mean")
	names(df_all_features)[ncol(df_all_features)]=paste0(feature,"_rank")

}
write.table(df_all_features, file="CNN_1000.mean_and_rank.per_feature.tab",sep="\t",quote=F)

## calculate p-values from normal distribution:

for ( i in 1:30){
	name=gsub("_mean","",names(df_all_features[2*i]))
	z=scale(df_all_features[,2*i], center=T,scale=T)
	p=2*(pnorm(-abs(z)))
	q=p.adjust(p,method="fdr")
	df_all_features=cbind(df_all_features, p)
	df_all_features=cbind(df_all_features, q)

	names(df_all_features)[ncol(df_all_features)-1]=paste0(name,"_p")
	names(df_all_features)[ncol(df_all_features)]=paste0(name,"_q")

}

colSums(df_all_features[,grep("_p",names(df_all_features))]<0.05)
# open_chromatin_p             tf_p       promoter_p       enhancer_p      repressed_p         active_p
#             4233             2210             4925             5195             5071             4019

colSums(df_all_features[,grep("_q",names(df_all_features))]<0.05)
# open_chromatin_q             tf_q       promoter_q       enhancer_q      repressed_q         active_q
#             1221              952             1735             1486             1455             1072

### how many variants significant for at least one group:            
length(which(rowSums(df_all_features[,grep("_q",names(df_all_features))]<0.05)>=1))
[1] 11389
length(which(rowSums(df_all_features[,grep("_q",names(df_all_features))]<0.01)>=1))
[1] 8850

### keep lowest Q per variant
df_all_features$lowest_Q=apply(df_all_features[,grep("_q",names(df_all_features))],1,min)
write.table(df_all_features, file="CNN_1000.mean_and_rank.p_and_q.per_feature.tab",sep="\t",quote=F)

### and lowest Q per feature group:
data$lowest_open_Q=apply(data[,c(62,64,66,72,80,88)],1,min)
data$lowest_TF_Q=apply(data[,c(78,86,112, 114, 116,118)],1,min)
data$lowest_promoter_Q=apply(data[,c(70,76,84,102,108,120)],1,min)
data$lowest_enhancer_Q=apply(data[,c(90, 94, 100, 110)],1,min)
data$lowest_repressed_Q=apply(data[,c(68,74,82,96,106)],1,min)
data$lowest_active_Q=apply(data[,c(92,98,104)],1,min)
write.table(data, file="CNN_1000.mean_and_rank.p_and_q.per_feature.tab",sep="\t",quote=F)


### how many variants are significant if we require Q<0.05 and absolute score difference >=0.05
colnames(data)[129] = "ASE_Q"
mean_cols = grep("mean",colnames(data))
q_cols = grep("_q",colnames(data))

select=unique(unlist(sapply(1:30, function(x) which(abs(data[,mean_cols[x]])>=0.05 & data[,q_cols[x]]<0.05))))
length(select)
# [1] 3464

data_signif=data[select,]
# set all the q's to 0.1 if score diff <0.05

for ( i in 1:30){
	data_signif[which(abs(data_signif[,mean_cols[i]])<0.05 & data_signif[,q_cols[i]]<0.05),q_cols[i]] = NA
}
### and lowest Q per feature group:
data_signif$lowest_Q=apply(data_signif[,grep("_q",names(data_signif))],1, function(x) min(x, na.rm=T))

data_signif$lowest_open_Q=apply(data_signif[,c(62,64,66,72,80,88)],1,function(x) min(x, na.rm=T))
data_signif$lowest_TF_Q=apply(data_signif[,c(78,86,112, 114, 116,118)],1,function(x) min(x, na.rm=T))
data_signif$lowest_promoter_Q=apply(data_signif[,c(70,76,84,102,108,120)],1,function(x) min(x, na.rm=T))
data_signif$lowest_enhancer_Q=apply(data_signif[,c(90, 94, 100, 110)],1,function(x) min(x, na.rm=T))
data_signif$lowest_repressed_Q=apply(data_signif[,c(68,74,82,96,106)],1,function(x) min(x, na.rm=T))
data_signif$lowest_active_Q=apply(data_signif[,c(92,98,104)],1,function(x) min(x, na.rm=T))



### barplot - ugly
bar=data.frame(colSums(data_signif[,c(121,138:143)]<0.05))
rownames(bar) = c("Any feature","Open chromatin","TF","Promoter","Enhancer","Repressed","Active")
names(bar)="Variants"
bar$Group=rownames(bar)
library(ggplot2)
pdf("signif_vars_by_category.pdf")
ggplot(data=bar, aes(x=Group, y=Variants)) +
  geom_bar(stat="identity") + scale_fill_brewer(palette = "Set1") + theme_bw()
dev.off()

#### make a barplot with 30bars, grouped into 6 categories
q_cols = grep("_q",colnames(data))
bar=data.frame(colSums(data[,q_cols]<0.05))

rownames(bar) = gsub("_q","",rownames(bar))
names(bar)="Variants"
bar$Feature=rownames(bar)
bar$Group=c(rep("Open chromatin",3), "Repressed", "Promoter","Open chromatin",  "Repressed", "Promoter", "TF","Open chromatin","Repressed", "Promoter","TF",
"Open chromatin", "Enhancer", "Active", "Enhancer", "Repressed", "Active", "Enhancer", "Promoter", "Active", "Repressed", "Promoter", "Enhancer", rep("TF",4), "Promoter" )
bar=bar[,c(3,2,1)]
bar$Group=factor(bar$Group, levels=c("Promoter","Enhancer", "Open chromatin", "Active", "TF", "Repressed"))
bar=bar[order(bar$Group),]
library(ggplot2)

pdf("signif_vars_by_category.pdf")
ggplot(bar, aes(x=Feature, y=Variants, fill=Group)) + 
  geom_bar( stat="identity", position=position_dodge(width=0.8)) +    
  scale_fill_manual(values=c("red","orange","yellow","lightgreen","cornflowerblue","grey")) + theme_bw(base_size = 14) + 
   xlab("Pancreatic islet ChromHMM state") + ylab("-log10(Q-value)")+ guides(fill = guide_legend(reverse=TRUE)) + coord_flip()

group_cols=c("red","orange","yellow","lightgreen","cornflowerblue","grey")
sapply(bar$Group, function(x) group_cols[which(levels(bar$Group)==x)])
pdf("signif_vars_by_category.pdf", width=10, height=5)
par(las=2, mar=c(10,5.1, 1,1))
barplot(bar$Variants, col=sapply(bar$Group, function(x) group_cols[which(levels(bar$Group)==x)]), border=NA, names.arg=bar$Feature, ylab="Variants with Q<0.05")
dev.off()

ggplot(chrom_plot, aes(x=variable, y=value, fill=Feature.group)) + 
  geom_bar( stat="identity", position=position_dodge(width=0.8)) +    
  scale_fill_manual(values=rev(c("red","orange","yellow","lightgreen","cornflowerblue","grey"))) + theme_bw(base_size = 14) + 
   xlab("Pancreatic islet ChromHMM state") + ylab("-log10(Q-value)")+ guides(fill = guide_legend(reverse=TRUE)) + coord_flip()



pdf("CNN_score_distribution.pdf", height=5, width=5)
par(mar=c(6,6,3,1))
for ( i in 1:length(mean_cols)){
	plot(density(data[,mean_cols[i]]), main=colnames(data[mean_cols[i]]), lwd=2, cex.axis=2, cex.lab=2,xlab="CNN score difference")
	abline(v=head(sort(data[which(data[,q_cols[i]]>=0.05),mean_cols[i]]))[1], col="red", lty=2, lwd=2)
	abline(v=tail(sort(data[which(data[,q_cols[i]]>=0.05),mean_cols[i]]),n=1), col="red", lty=2, lwd=2)

}

dev.off()

