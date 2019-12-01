setwd("~/Desktop/islets CNN/GWAS_CNN_convergence/")

load("cred_sets.Rdata")

locus_ann=read.table("~/Desktop/islet T2D loci/HRC_locus_annotation.txt",h=T,sep="\t")
locus_ann$locus=paste(locus_ann$Nearest_gene,locus_ann$Index_variant,locus_ann$Novel, locus_ann$Index,sep="_")
locus_ann$IndexSNP=paste(locus_ann$Chr, locus_ann$Position, sep="_")
cred_names=do.call(rbind,lapply(cred, function(x)   subset(locus_ann, IndexSNP==x[1,]$IndexSNP, c(locus,Jason_signal_id))))

for (i in 1:length(cred)){
  credset=cred[[i]]
  fPPA0.2=nrow(subset(credset, FGWAS_PPA>=0.2))
}

## how many signals do we have with >=2 variants with fPPA>=0.2
length(which(unlist(lapply(cred, function(x) nrow(subset(x, FGWAS_PPA>=0.2))))>=2))
# 95

cred_fPPA0.2=lapply(cred, function(x) subset(x, FGWAS_PPA>=0.2))

### find how many signals are already fine-mapped to single variant by gPPA or fPPA (>=0.8)
### how many of these are predicted regulatory by CNNs

cred_PPA0.8=lapply(cred, function(x) subset(x, PPAg >=0.8 | FGWAS_PPA>=0.8))
sum(unlist(lapply(cred_PPA0.8, function(x) ifelse(nrow(x)>=1,1,0))))  ## 74
sum(unlist(lapply(cred_PPA0.8, function(x) ifelse((x$PPAg >=0.8 | x$FGWAS_PPA>=0.8) & x$lowest_Q<0.05,1,0))))  ## 15

cred_PPA0.8_df=do.call(rbind, cred_PPA0.8)
### Sup. Table 4- true positives
subset(cred_PPA0.8_df, lowest_Q<0.05)

sum(unlist(lapply(cred_PPA0.8, function(x) ifelse(x$PPAg >=0.8 ,1,0))))  ## 52
sum(unlist(lapply(cred_PPA0.8, function(x) ifelse(x$PPAg >=0.8 & x$lowest_Q<0.05,1,0))))  ## 15

sum(unlist(lapply(cred_PPA0.8, function(x) ifelse(x$FGWAS_PPA >=0.8 ,1,0))))  ## 69
sum(unlist(lapply(cred_PPA0.8, function(x) ifelse(x$FGWAS_PPA >=0.8 & x$lowest_Q<0.05,1,0))))  ## 28

### find signals with no variant fine-mapped above PPA>=0.8, but with multiple variants PPA>=0.2
### how many of these are predicted regulatory by CNNs

finemap_plot_input=list()
fine_mappable_signals=0
for (i in 1:length(cred)){
  if (nrow(subset(cred[[i]], PPAg>=0.8 | FGWAS_PPA>=0.8))==0){
    if (nrow(subset(cred[[i]], FGWAS_PPA>=0.2))>=2){
      fine_mappable_signals=fine_mappable_signals+1
      if (nrow(subset(cred[[i]], FGWAS_PPA>=0.2 & lowest_Q<0.05))>=1){
        print(paste0(cred_names[i,]$locus,"\t",nrow(subset(cred[[i]], FGWAS_PPA>=0.2)),"\t",
            nrow(subset(cred[[i]], FGWAS_PPA>=0.2 & lowest_Q<0.05))))
        finemap_plot_input[[i]]=data.frame(locus=rep(cred_names[i,]$locus, nrow(cred[[i]])), variant=cred[[i]]$variant, pos=cred[[i]]$Pos,
                      lowest_Q=cred[[i]]$lowest_Q, Genetic_PPA=cred[[i]]$PPAg, FGWAS_PPA=cred[[i]]$FGWAS_PPA,
                      ChromHMM=cred[[i]]$ChromHMM)
      }
    }
  }
}

fine_mappable_signals
## 93

### plots
finemap_plotting=do.call(rbind, finemap_plot_input)

length(levels(finemap_plotting$locus))
## 37

finemap_plotting$lowest_Q[finemap_plotting$lowest_Q==0] = 1e-250
finemap_plotting=data.frame(finemap_plotting)

finemap_plotting$locus_name_plot=sapply(strsplit(as.character(finemap_plotting$locus), split="_"), function(x) 
  paste(x[1],x[2]))

finemap_plotting$top_feature=rep(NA, nrow(finemap_plotting))
for(v in as.character(finemap_plotting$variant)){
 tmp=subset(data, variant==v & lowest_Q<0.05)
  if(nrow(tmp)>0){
    names=names(tmp)[which(tmp==tmp$lowest_Q)]
    names=names[grepl("_q",names)][1]
    finemap_plotting$top_feature[which(finemap_plotting$variant==v)]=gsub("_q","",names)
 }
}
  
finemap_plotting$location=sapply(as.character(finemap_plotting$variant), function(x) subset(data, variant==x, pos)[,1])
                                    
                                    
pdf("fine_map_plots.301119.pdf",width=10)
for (l in unique(finemap_plotting$locus_name_plot)){
  tmp=finemap_plotting[finemap_plotting$locus_name_plot==l & (finemap_plotting$Genetic_PPA>=0.01 | finemap_plotting$FGWAS_PPA>=0.01),]
  tmp=tmp[order(tmp$pos),]
  width=nrow(tmp)/2
  if (nrow(tmp)<10){width=5}
  par(mar=c(1.2,7,1,1), oma=c(8,1,4,1), mfrow=c(3,1), xpd=F, mgp=c(3.5,1,0))
  plot(1:length(tmp$variant), tmp$Genetic_PPA, pch=21,cex=2, col="black",bg="cornflowerblue", xaxt="n",xlab="",ylab="Genetic PPA", 
       cex.lab=1.5, cex.axis=1.5, ylim=c(0,max(tmp$Genetic_PPA)+0.1))
  plot(1:length(tmp$variant), tmp$FGWAS_PPA, pch=21, cex=2,col="black",bg="darkolivegreen3", xaxt="n",xlab="",ylab="Functional PPA", 
       cex.lab=1.5, cex.axis=1.5, ylim=c(0,max(tmp$FGWAS_PPA,na.rm=T)+0.1))
  plot(1:length(tmp$variant), -log10(tmp$lowest_Q), pch=21,cex=2, col="black",bg="coral", xaxt="n",xlab="",ylab="-log10 CNN lowest Q", 
       cex.lab=1.5, cex.axis=1.5, ylim=c(0,max(-log10(tmp$lowest_Q))+0.25))
  abline(h=-log10(0.05), lty=2, col="red")
  par(xpd=NA,las=2)
  axis(1,at=c(1:length(tmp$variant)),labels=paste0(tmp$variant,"\n",tmp$location), tick=F, cex.axis=1.5)
  ## add text with best feature for the min variant:
  text(x=which.min(tmp$lowest_Q), y=-log10(tmp$lowest_Q[which.min(tmp$lowest_Q)])*0.9, 
  tmp$top_feature[which.min(tmp$lowest_Q)])
  
  par(las=1)
  mtext(l, outer = TRUE, cex = 1.5)
}
dev.off()
                                 
### Sup. Table 5
subset(finemap_plotting, FGWAS_PPA>=0.05)

### get the top CNN feature and the mean difference for it
data=data.frame(fread("CNN_1000.mean_and_rank.p_and_q.per_feature.tab")) ### this is variant level
colnames(data)[122] = "lowest_Q"
colnames(data)[123] = "Genetic_PPA"
names(data)[130]="ASEQ_Q"

vars=as.character(subset(finemap_plotting, FGWAS_PPA>=0.05,variant)[,1])
top_CNN_res=data.frame(variant=vars, top_CNN_feature=rep(NA,length(vars)), CNN_score_diff=rep(0,length(vars)))

for (i in 1:length(vars)){
  r=rownames(subset(data, V1==vars[i]))
  top_CNN_res$top_CNN_feature[i]=gsub("_q","",names(sort(data[r,grepl("_q", colnames(data))])[1]))
  top_CNN_res$CNN_score_diff[i]=data[r,paste0(top_CNN_res$top_CNN_feature[i],"_mean")]
}
