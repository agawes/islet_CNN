## read in my results && get E-values
setwd("/well/got2d/agata/Basset/Full_islet_dataset")
data=read.table("CNN_1000.mean_and_rank.p_and_q.per_feature.tab")
options(width=200)


## read in all DeepSEA results
setwd("/well/got2d/agata/Basset/HRC_credset/DeepSEA")
deepsea_files=list.files(".",pattern="infile.vcf.out.funsig", recursive=T)

deepsea_all_res=data.frame()
for (I in deepsea_files){
	tmp=read.csv(I,h=T)
	deepsea_all_res=rbind(deepsea_all_res,tmp)
}
deepsea_all_res=deepsea_all_res[which(!duplicated(deepsea_all_res$name)),]
rownames(deepsea_all_res)=deepsea_all_res$name

merged=merge(data, deepsea_all_res,by="row.names")
cor(merged$lowest_Q, merged$DeepSEA.score)
# [1] 0.2270982
cor.test(merged$lowest_Q, merged$DeepSEA.score)$p.value
# p-value < 2.2e-16

merged$CNN_signif=ifelse(merged$lowest_Q<0.05,1,0)
merged$DeepSEA_signif=ifelse(merged$DeepSEA.score<0.05,1,0)

table(merged[,c(152,153)])
          DeepSEA_signif
CNN_signif     0     1
         0 85872 12041
         1  7304  4068
         

merged$col=rep(rgb(0,0,0,50,maxColorValue=255),nrow(merged))
merged$col[merged$lowest_Q<0.05 & merged$DeepSEA.scor<0.05] = rgb(100,0,0,50,maxColorValue=255)
merged$col[merged$lowest_Q<0.05 & merged$DeepSEA.scor>=0.05] = rgb(0,100,0,50,maxColorValue=255)
merged$col[merged$lowest_Q>=0.05 & merged$DeepSEA.scor<0.05] = rgb(0,0,100,50,maxColorValue=255)


pdf("CNN_vs_DeepSEA.pdf")
par(mar=c(5.1,5.1,1,1))
plot(-log10(merged$lowest_Q), -log10(merged$DeepSEA.score), col=merged$col, pch=16, bty="n", xlab="-log10(lowest Q) islet CNNs", ylab="-log10(funsig) DeepSEA", cex.axis=1.6, cex.lab=1.6 )
text(25,5, "r=0.227", cex=2)
dev.off()


deepsea_panc_files=list.files(".",pattern="infile.vcf.out.evalue", recursive=T)

deepsea_panc_res=data.frame()
for (I in deepsea_panc_files){
	#tmp=read.csv(I,h=T)[,c(1:6,38)]
	tmp=read.csv(I,h=T,colClasses=c(rep(NA,6),rep("NULL",31),NA,rep("NULL",887)) )
	deepsea_panc_res=rbind(deepsea_panc_res,tmp)
}
deepsea_panc_res=deepsea_panc_res[which(!duplicated(deepsea_panc_res$name)),]
rownames(deepsea_panc_res)=deepsea_panc_res$name

merged_panc=merge(data, deepsea_panc_res,by="row.names")
cor(merged_panc$lowest_Q, merged_panc$PanIslets.DNase.None)
[1] 0.223
cor(-log10(merged_panc$lowest_Q), -log10(merged_panc$PanIslets.DNase.None))

merged_panc$CNN_signif=ifelse(merged_panc$lowest_Q<0.05,1,0)
merged_panc$DeepSEA_signif=ifelse(merged_panc$PanIslets.DNase.None<0.05,1,0)


table(merged_panc[,c(152,153)])
          DeepSEA_signif
CNN_signif     0     1
         0 92811  5102
         1  8387  2985
                  
merged_panc$col=rep(rgb(0,0,0,50,maxColorValue=255),nrow(merged_panc))
merged_panc$col[merged_panc$lowest_Q<0.05 & merged_panc$PanIslets<0.05] = rgb(100,0,0,50,maxColorValue=255)
merged_panc$col[merged_panc$lowest_Q<0.05 & merged_panc$PanIslets>=0.05] = rgb(0,100,0,50,maxColorValue=255)
merged_panc$col[merged_panc$lowest_Q>=0.05 & merged_panc$PanIslets<0.05] = rgb(0,0,100,50,maxColorValue=255)


pdf("CNN_vs_DeepSEA.pdf",width=14)
par(mar=c(5.1,5.1,1,1),mfrow=c(1,2))
plot(-log10(merged$lowest_Q), -log10(merged$DeepSEA.score), col=merged$col, pch=16, bty="n", xlab="-log10(lowest Q) islet CNNs", ylab="-log10(funsig) DeepSEA", cex.axis=1.6, cex.lab=1.6, ylim=c(0,5.5) )
text(25,5, "r=0.227", cex=2)

plot(-log10(merged_panc$lowest_Q), -log10(merged_panc$PanIslets.DNase.None), col=merged_panc$col, pch=16, bty="n", xlab="-log10(lowest Q) islet CNNs", ylab="-log10(E-value) DeepSEA: PanIslets DNase", cex.axis=1.6, cex.lab=1.6, ylim=c(0,5.5) )
text(25,5, "r=0.223", cex=2)
dev.off()

    

