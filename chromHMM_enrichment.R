## chrom HMM enrichment plots
library(data.table)
library(gage)
library(ggplot2)

setwd("~/Desktop/islets CNN/GWAS_CNN_convergence/")

data=data.frame(fread("CNN_1000.mean_and_rank.p_and_q.per_feature.tab"))
colnames(data)[122] = "lowest_Q"
colnames(data)[123] = "Genetic_PPA"

data$ChromHMM=factor(data$ChromHMM)
snps_by_state = list()
## for all chromHMM states:
for(i in levels(data$ChromHMM)) {
  snps_by_state[[i]] = as.character(rownames(data[which(data$ChromHMM==i),]))
}


chrom_df=data.frame(matrix(vector(), 0, length(levels(data$ChromHMM)), dimnames=list(c(),levels(data$ChromHMM))))
indexes=grep("lowest",names(data))
for (i in 1:length(indexes)){
  input=data[,indexes[i]]
  names(input) = rownames(data)
  gage_res <- gage(input, gsets = snps_by_state, ref = NULL, samp = NULL, same.dir = T, set.size = c(10, 41000), full.table=T)
  row=gage_res$less[levels(data$ChromHMM),4]
  chrom_df[i,]=row
}
rownames(chrom_df)=names(data)[grep("lowest",names(data))]

### keep the larger of enrichments
chrom_df=data.frame(matrix(vector(), 0, length(levels(data$ChromHMM)), dimnames=list(c(),levels(data$ChromHMM))))
indexes=grep("lowest",names(data))
for (i in 1:length(indexes)){
  input=data[,indexes[i]]
  names(input) = rownames(data)
  gage_res <- gage(input, gsets = snps_by_state, ref = NULL, samp = NULL, same.dir = T, set.size = c(10, 41000), full.table=T)
  gage_top=rep(0,length(levels(data$ChromHMM)))
  names(gage_top)=levels(data$ChromHMM)
  for (l in levels(data$ChromHMM)){
    gage_top[l]=ifelse(gage_res$less[l,4]<gage_res$greater[l,4],-log10(gage_res$less[l,4]), log10(gage_res$greater[l,4]))
  }
  chrom_df[i,]=gage_top
}
rownames(chrom_df)=names(data)[grep("lowest",names(data))]

library(reshape)
chrom_df$Group=rownames(chrom_df)
chrom_melt=melt(chrom_df)

chrom_plot=chrom_melt[with(chrom_melt, variable %in% c("E11","E12","E8","E9","E4","E5","E2") & Group != "lowest_Q" & Group != "log10_lowestQ_plot"),]
chrom_plot$variable=factor(chrom_plot$variable, levels=c("E2","E5","E4","E9","E8","E12","E11"))
levels(chrom_plot$variable) = c("Low Methylation", "Heterochromatin", "Lowly-methylated\nWeak Enhancer","Open Weak Enhancer","Open Strong Enhancer","Weak Promoter","Active Promoter")
chrom_plot$Group=factor(chrom_plot$Group, levels=c("lowest_repressed_Q","lowest_TF_Q","lowest_open_Q", "lowest_enhancer_Q","lowest_active_Q","lowest_promoter_Q"))
levels(chrom_plot$Group) = c("Repressed", "TF","Open chromatin","Enhancer","Active","Promoter")
chrom_plot$Group=factor(chrom_plot$Group, levels=c("Repressed", "TF","Active","Open chromatin","Enhancer","Promoter"))

names(chrom_plot)[1] = "Feature.group"
# pdf("chromHMM_enrichments.bar.v.pdf") # ugly
# ggplot(chrom_plot, aes(x=variable, y=value, fill=Feature.group)) + 
#   geom_bar( stat="identity", position=position_dodge(width=0.8)) + coord_flip() +
#   scale_fill_manual(values=rev(c("red","orange","yellow","lightgreen","cornflowerblue","grey"))) + theme_bw(base_size = 14) + 
#   xlab("Pancreatic islet ChromHMM state") + ylab("-log10(Q-value)")+ guides(fill = guide_legend(reverse=TRUE)) 
# dev.off()

chrom_plot$Feature.group=factor(chrom_plot$Feature.group, levels=rev(c("Repressed", "TF","Active","Open chromatin","Enhancer","Promoter")))
chrom_plot$variable=factor(chrom_plot$variable, levels=rev(c("Low Methylation", "Heterochromatin", "Lowly-methylated\nWeak Enhancer","Open Weak Enhancer","Open Strong Enhancer","Weak Promoter","Active Promoter")))

text_format <- element_text(face = "bold", color="black",size = 16)

pdf("chromHMM_enrichments.bar.h.pdf", width=12, height=6) # ugly
ggplot(chrom_plot, aes(x=variable, y=value, fill=Feature.group)) + 
  geom_bar( stat="identity", position=position_dodge(width=0.8)) + 
  scale_fill_manual(values=c("red","orange","yellow","lightgreen","cornflowerblue","grey")) + theme_bw(base_size = 14) + 
  xlab("Pancreatic islet ChromHMM state") + ylab("-log10(Q-value)")+ guides(fill = guide_legend()) +
  theme(axis.text=text_format,axis.title=text_format, legend.text=text_format,  legend.title = text_format,
        axis.text.x = element_text(angle = 30, hjust = 1), legend.position = c(0.75, 0.75))
dev.off()


####  5/07/19 - make a heatmap instead of barplots:
colnames(chrom_df)=c("PR","OC","AP","WP","GE","TE","INS","LMR","HET2","LMWE","HET","CWE","CSE","OSE","OWE")
library(reshape)
chrom_df$Group=rownames(chrom_df)
chrom_melt=melt(chrom_df)

chrom_plot=chrom_melt[with(chrom_melt, Group != "lowest_Q" & Group != "log10_lowestQ_plot"),]
chrom_plot$Group=factor(chrom_plot$Group, levels=c("lowest_repressed_Q","lowest_TF_Q","lowest_open_Q", "lowest_enhancer_Q","lowest_active_Q","lowest_promoter_Q"))
levels(chrom_plot$Group) = c("Repressed", "TF","Open chromatin","Enhancer","Active","Promoter")
chrom_plot$Group=factor(chrom_plot$Group, levels=c("Repressed", "TF","Active","Open chromatin","Enhancer","Promoter"))
chrom_plot$q_star=cut(chrom_plot$value, breaks=c(-Inf, -3, -2, -1, 1,2,3,Inf), 
                      label=c("***", "**", "*", "", "*","**","***")) 
chrom_plot$variable=factor(chrom_plot$variable, levels=c("AP","WP","OSE","OWE","CSE","CWE","LMWE","GE",
                                                         "OC","PR","INS","TE","HET2","HET","LMR"))

### trim all values at +/- 100 
chrom_plot$value[chrom_plot$value<=(-50)] = -50
chrom_plot$value[chrom_plot$value>50] = 50

pdf("ChromHMM.enrichment.heatmap.pdf", width=10, height=5)
text_format <- element_text( color="black",size = 16)
p <- ggplot(aes(x=variable, y=Group, fill=value), data=chrom_plot) + geom_tile() + 
  scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
 # geom_text(aes(label=q_star), color="black", size=3) + 
  labs(y="CNN feature group", x="Islet chromatin states", fill="Enrichment\n-log10(p)") + theme_bw() + 
  theme(axis.text.x=element_text(angle = 90, hjust = 0, color="black",size=16),
      axis.text.y=text_format, axis.title=text_format) 
  
p
dev.off()
