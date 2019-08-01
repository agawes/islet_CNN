library(data.table)
setwd("~/Desktop/islets CNN/GWAS_CNN_convergence/")

data=data.frame(fread("CNN_1000.mean_and_rank.p_and_q.per_feature.tab")) ### this is variant level
colnames(data)[122] = "lowest_Q"
colnames(data)[123] = "Genetic_PPA"

fgwas=read.table("~/Desktop/islet T2D loci/islet_fGWAS_PPAs.results_func-cred-sets.txt",h=T,sep="\t")
## ann
locus_ann=read.table("~/Desktop/islet T2D loci/HRC_locus_annotation.txt",h=T,sep="\t")
locus_ann$locus=paste(locus_ann$Nearest_gene,locus_ann$Index_variant,locus_ann$Novel, locus_ann$Index,sep="_")


#### I need to have this all stored as list, to preserve the credible sets structure:
## first read in all the results
credset_dir="HRC_per_locus_credsets"
files<-list.files(path=credset_dir, pattern=".txt$",recursive=TRUE)

cred=list()
for (f in files){
  name = gsub("credible_set_Eur_","",gsub(".txt","",f))
  df = read.table(paste0(credset_dir,"/",f),header=T,sep="\t")
  cred[[name]] = df
  print(f)
}

### get a matching signal id data frame
locus_ann$IndexSNP=paste(locus_ann$Chr, locus_ann$Position, sep="_")
cred_names=do.call(rbind,lapply(cred, function(x)   subset(locus_ann, IndexSNP==x[1,]$IndexSNP, c(locus,Jason_signal_id))))

### append lowest_Q, and FGWAS_PPA, ChromHMM, rsID
for (sign in 1:length(cred)){
  cred[[sign]]$variant=rep("",nrow(cred[[sign]]))
  cred[[sign]]$ChromHMM=rep(NA,nrow(cred[[sign]]))
  cred[[sign]]$lowest_Q=rep(1,nrow(cred[[sign]]))
  cred[[sign]]$FGWAS_PPA=rep(0,nrow(cred[[sign]]))

  for (var in 1:nrow(cred[[sign]])){
    position=paste0(cred[[sign]][var,]$Chr,":", cred[[sign]][var,]$Pos)
    
    cred[[sign]]$variant[var]=subset(data, pos==position, V1)[1,]
    chrom=subset(data, pos==position, ChromHMM)[1,]
    if (!is.na(chrom)){cred[[sign]]$ChromHMM[var]=as.character(chrom)}
    
    q=subset(data, pos==position, lowest_Q)[1,]
    if (!is.na(q)){
      cred[[sign]]$lowest_Q[var]=as.numeric(q)
    }
  fPPA=subset(fgwas, SNPID==paste0("chr",position) & Locus.ID==as.character(cred_names[sign,]$Jason_signal_id), PPA)
  if (nrow(fPPA)==1){
       cred[[sign]]$FGWAS_PPA[var]=as.numeric(fPPA)
  }
  if (nrow(fPPA)>1){
    print(paste(position,"problematic - >1 fPPA"))
  }
  }
  print(sign)
}

save(cred, file="cred_sets.Rdata")
