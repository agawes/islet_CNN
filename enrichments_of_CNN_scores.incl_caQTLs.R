#### analysis needed for Figure 1 in CNN paper
library(gage)
library(ggplot2)

data=read.table("CNN_1000.mean_and_rank.p_and_q.per_feature.tab",h=T)
variants=gsub("_alt","",data$variant)
rownames(data)=variants

data=data[,-1]

### add the PPAs, chromHMM, etc
ppa=read.table("rank_product.1000nets.all_SNPs.fGWAS_PPA.txt",h=T,sep="\t")
ppa=ppa[!duplicated(ppa$V1),]
rownames(ppa)=ppa$V1
names(ppa)[3:4] = c("Genetic_PPA","FGWAS_PPA")
data=cbind(data, ppa[variants,3:4])

### FGWAS_PPA annotation is wrong for some signals were SNPs may be present in the file twice
### let's double check this with the original file from Jason
fgwas_ann=read.table("../HRC_credset/islet_fGWAS_PPAs.results_func-cred-sets.txt",h=T,sep="\t")
fgwas_ppas=sapply(as.character(data$pos), function(x) subset(fgwas_ann, SNPID==paste0("chr",x),PPA))

problematic=which(sapply(fgwas_ppas, length)>=2)



chromHMM=read.table("../HRC_credset/HRC_SNPs.chromHMM_states.txt")
chromHMM=chromHMM[!duplicated(chromHMM$V1),]
rownames(chromHMM)=chromHMM$V1
names(chromHMM)[2] = "ChromHMM"
data=cbind(data, chromHMM[variants,2,drop=F])

ase=read.table("../HRC_credset/HRC.ASE_Elife2018.txt",h=T)
ase=ase[!duplicated(ase$snp),]
rownames(ase)=ase$snp
names(ase)=c("snp","ASE_al1","ASE_al2","ASE_ratio","ASE_p")
ase$ASE_q=p.adjust(ase$ASE_p, method="fdr")
data=cbind(data, ase[variants,2:6])

ann_gerp=read.table("rank_product.1000nets.all_SNPs.ann_gerp.txt",h=T,sep="\t")
ann_gerp=ann_gerp[!duplicated(ann_gerp$V1),]
rownames(ann_gerp)=ann_gerp$V1
names(ann_gerp)[10]="locus"
data=cbind(data, ann_gerp[variants,10,drop=F])

write.table(data, file="CNN_1000.mean_and_rank.p_and_q.ann.per_feature.tab",sep="\t",quote=F)



mean_cols = grep("mean",colnames(data))
q_cols = grep("_q",colnames(data))

### keep lowest Q per variant
colnames(data)[129] = "ASE_Q"
data$lowest_Q=apply(data[,grep("_q",names(data))],1,function(x) min(x, na.rm=T))

### and lowest Q per feature group:
data$lowest_open_Q=apply(data[,c(62,64,66,72,80,88)],1,function(x) min(x, na.rm=T))
data$lowest_TF_Q=apply(data[,c(78,86,112, 114, 116,118)],1,function(x) min(x, na.rm=T))
data$lowest_promoter_Q=apply(data[,c(70,76,84,102,108,120)],1,function(x) min(x, na.rm=T))
data$lowest_enhancer_Q=apply(data[,c(90, 94, 100, 110)],1,function(x) min(x, na.rm=T))
data$lowest_repressed_Q=apply(data[,c(68,74,82,96,106)],1,function(x) min(x, na.rm=T))
data$lowest_active_Q=apply(data[,c(92,98,104)],1,function(x) min(x, na.rm=T))


#############################################################################



################  check enrichment of caQTLs in the 6 ranked lists: ########################


ase_snps = list()
ase_snps[["p0.05"]] = as.character(ase$snp[which(ase$ASE_p<0.05)])
ase_snps[["q0.05"]]= as.character(ase$snp[which(ase$ASE_q<0.05)])

ase_df=setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("ASE_p", "ASE_q"))

for (i in c(121, 138:143)){
	input=data[,i]
	names(input) = rownames(data)

	gage_res <- gage(input, gsets = ase_snps, ref = NULL, samp = NULL, same.dir = T, set.size = c(10, 41000), full.table=T)
	ase_df=rbind(ase_df,gage_res$less[,4])
}
rownames(ase_df)=names(data)[c(121, 138:143)]
names(ase_df)=c("ASE_p", "ASE_q")

### 
#                 ASE_p        ASE_q
# lowest_Q 6.695492e-33 6.925428e-20



################  enrichment of chromHMM states ###########

snps_by_state = list()
## for all chromHMM states:
for(i in levels(data$ChromHMM)) {
    snps_by_state[[i]] = as.character(rownames(data[which(data$ChromHMM==i),]))
}

chrom_df=setNames(data.frame(matrix(ncol = length(levels(data$ChromHMM)), nrow = 0)), levels(data$ChromHMM))

#for ( i in 61:121){
for (i in c(121, 138:143)){

	input=data[,i]
	names(input) = rownames(data)
	gage_res <- gage(input, gsets = snps_by_state, ref = NULL, samp = NULL, same.dir = T, set.size = c(10, 41000), full.table=T)
	chrom_df=rbind(chrom_df,gage_res$less[levels(data$ChromHMM),4])
}
#rownames(chrom_df)=names(data)[61:121]
rownames(chrom_df)=names(data)[c(121, 138:143)]

names(chrom_df)=levels(data$ChromHMM)

 chrom_df[61,]
#         E1       E10           E11          E12       E13 E14         E15         E2 E3           E4 E5        E6 E7           E8           E9
#  lowest_Q  1 0.4907798 6.281951e-136 1.102186e-31 0.9879558   1 0.001193512 0.09253855  1 1.702104e-08  1 0.8422992  1 2.305458e-51 2.310253e-25

### keep the larger of enrichments
chrom_df=setNames(data.frame(matrix(ncol = length(levels(data$ChromHMM)), nrow = 0)), levels(data$ChromHMM))
for (i in c(120, 137:142)){
	input=data[,i]
	names(input) = rownames(data)
	gage_res <- gage(input, gsets = snps_by_state, ref = NULL, samp = NULL, same.dir = T, set.size = c(10, 41000), full.table=T)
	gage_top=rep(0,length(levels(data$ChromHMM)))
	names(gage_top)=levels(data$ChromHMM)
	for (l in levels(data$ChromHMM)){
		gage_top[l]=ifelse(gage_res$less[l,4]<gage_res$greater[l,4],-log10(gage_res$less[l,4]), log10(gage_res$greater[l,4]))
	}
	chrom_df=rbind(chrom_df,gage_top)
}
rownames(chrom_df)=names(data)[c(120, 137:142)]
names(chrom_df)=levels(data$ChromHMM)

#pdf("chromHMM_enrichments.heatmap.pdf") # ugly
#heatmap(data.matrix(chrom_df[,1:15]),col=brewer.pal(9,"RdBu"))

#ggplot(chrom_melt, aes(x=variable, y=Group, fill=value)) + 
                  geom_tile()
#dev.off()


library(reshape)
chrom_df$Group=rownames(chrom_df)
chrom_melt=melt(chrom_df)

chrom_plot=chrom_melt[with(chrom_melt, variable %in% c("E11","E12","E8","E9","E4","E5","E2") & Group != "lowest_Q"),]
chrom_plot$variable=factor(chrom_plot$variable, levels=c("E2","E5","E4","E9","E8","E12","E11"))
levels(chrom_plot$variable) = c("Low Methylation", "Heterochromatin", "Lowly-methylated\nWeak Enhancer","Open Weak Enhancer","Open Strong Enhancer","Weak Promoter","Active Promoter")
chrom_plot$Group=factor(chrom_plot$Group, levels=c("lowest_repressed_Q","lowest_TF_Q","lowest_open_Q", "lowest_enhancer_Q","lowest_active_Q","lowest_promoter_Q"))
levels(chrom_plot$Group) = c("Repressed", "TF","Open chromatin","Enhancer","Active","Promoter")
chrom_plot$Group=factor(chrom_plot$Group, levels=c("Repressed", "TF","Active","Open chromatin","Enhancer","Promoter"))

names(chrom_plot)[1] = "Feature.group"
pdf("chromHMM_enrichments.bar.v.pdf") # ugly
ggplot(chrom_plot, aes(x=variable, y=value, fill=Feature.group)) + 
  geom_bar( stat="identity", position=position_dodge(width=0.8)) + coord_flip() +
  scale_fill_manual(values=rev(c("red","orange","yellow","lightgreen","cornflowerblue","grey"))) + theme_bw(base_size = 14) + 
   xlab("Pancreatic islet ChromHMM state") + ylab("-log10(Q-value)")+ guides(fill = guide_legend(reverse=TRUE)) 
dev.off()

chrom_plot$Feature.group=factor(chrom_plot$Feature.group, levels=rev(c("Repressed", "TF","Active","Open chromatin","Enhancer","Promoter")))

pdf("chromHMM_enrichments.bar.h.pdf", width=12, height=6) # ugly
ggplot(chrom_plot, aes(x=variable, y=value, fill=Feature.group)) + 
  geom_bar( stat="identity", position=position_dodge(width=0.8)) + 
  scale_fill_manual(values=c("red","orange","yellow","lightgreen","cornflowerblue","grey")) + theme_bw(base_size = 14) + 
   xlab("Pancreatic islet ChromHMM state") + ylab("-log10(Q-value)")+ guides(fill = guide_legend()) 
dev.off()


chrom_melt$variable=factor(chrom_melt$variable, levels=c("E2","E5","E3","E1","E15","E14","E13","E10","E6","E7","E4","E9","E8","E12","E11"))
chrom_melt$Group=factor(chrom_melt$Group, levels=c("lowest_repressed_Q","lowest_TF_Q","lowest_open_Q", "lowest_enhancer_Q","lowest_active_Q","lowest_promoter_Q"))
levels(chrom_melt$Group) = c("Repressed", "TF","Open chromatin","Enhancer","Active","Promoter")

#pdf("chromHMM_enrichments.all.bar.1.pdf") # ugly
#ggplot(chrom_melt, aes(x=variable, y=value, fill=Group)) + 
#  geom_bar( stat="identity", position=position_dodge(width=0.8)) +    
#  scale_fill_brewer(palette = "Set1") + theme_bw(base_size = 14) + 
#   xlab("") + guides(fill = guide_legend(reverse=TRUE)) + coord_flip()
#dev.off()



### twister plot  - OLD

chrom=read.table("chromHMM.tornado_plot.inout",h=T,sep="\t")
chrom$col=rep("grey",nrow(chrom))
chrom$col[which(chrom$X.log10_P>(-log10(0.05)))]="red"
chrom$col[which(chrom$X.log10_P<(log10(0.05)))]="blue"

chrom$hjust <- ifelse(chrom$X.log10_P > 0, -30, 26)
pdf("twister.chromHMM_states.pdf",width=15)
par(las=2, mar=c(5,5,3,1),xpd=NA)
bar=barplot(chrom$X.log10_P, ylab="-log10(P-value)", col=chrom$col, cex.axis=1.5, cex.lab=1.5)

text(bar,chrom$hjust, chrom$Annotation, srt=90, cex=1.4)
dev.off()





