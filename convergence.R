library(data.table)
library(nnet) ## which.is.max function
library(gage)
library("DescTools")

setwd("~/Desktop/islets CNN/GWAS_CNN_convergence/")

data=data.frame(fread("CNN_1000.mean_and_rank.p_and_q.per_feature.tab"))
colnames(data)[122] = "lowest_Q"
colnames(data)[123] = "Genetic_PPA"


load("cred_sets.Rdata")

### how often is the top variant by PPA significant with CNNs?
sum(unlist(lapply(cred, function(x) ifelse(x[1,]$lowest_Q<0.05,1,0))))  ### 61  (16%)

### we could get the random chance of getting this result by permuting the CNN p-values, but keeping the locus structure
### let's do this 1000 times, and take the mean, or plot all distribution?
perm_data=subset(data, , c(pos, lowest_Q, Genetic_PPA, FGWAS_PPA))
n_perm=1000
set.seed(1234)
for (i in 1:n_perm){
  perm_data[,4+i]=sample(perm_data$lowest_Q)
}

 
cred_perm=cred
for (c in 1:length(cred)){
    print(c)
  
    ### append the CNN permuted Q-values
    tmp=matrix(,nrow=nrow(cred[[c]]), ncol=n_perm)  
    for (i in 1:nrow(cred[[c]])){
        loc=paste(cred[[c]][i,]$Chr, cred[[c]][i,]$Pos,sep=":")
        tmp[i,]=as.numeric(subset(perm_data, pos==loc )[1,c(3:(n_perm+2))])
    }
    cred_perm[[c]]=cbind(cred[[c]], tmp)
}

library(ggplot2)
text_format <- element_text(face = "bold", color="black",size = 16)
cols <- c("Original"="black","Permuted"="steelblue")

#### going down in PPA from 1.00 to 0.01 in steps of 0.01 - what % of variants are significant <0.05 at this or higher PPA
cred_df=rbindlist(cred)
signif_fraction=sapply(seq(from=1,to=0.00,by=-0.01), function(x) 
  length(which(subset(cred_df, PPAg>=x, lowest_Q)[,1]<0.05))/nrow(subset(cred_df, PPAg>=x)))
plot(1:101,signif_fraction, ylim=c(0,0.6), type="l")

AUC_gPPA=AUC(1:101, signif_fraction)

### add line for random:
cred_perm_df=rbindlist(cred_perm)
cred_perm_df=data.frame(cred_perm_df)
#cred_perm_df=cred_perm_df[!(duplicated(paste0(cred_perm_df$Chr,":",cred_perm_df$Pos))),]

signif_fraction_perm=matrix(,nrow=101,ncol=n_perm)
for (i in 1:n_perm){
  signif_fraction_perm[,i]=sapply(seq(from=1,to=0.00,by=-0.01), function(x) 
    length(which(subset(cred_perm_df, PPAg>=x,i+5)[,1]<0.05))/nrow(subset(cred_perm_df, PPAg>=x)))
  print(i)
}

## p-value of this enrichment:
length(which(apply(signif_fraction_perm, 2, function(x) AUC(1:101,x)) >= AUC_gPPA))/n_perm

### get the enrichment p-value - see in how many permutations AUC is higher than of the real CNN q-values

ppa_df=data.frame(PPA=seq(from=1,to=0.00,by=-0.01), original=signif_fraction,
                  perm_mean=sapply(1:101, function(x) mean(signif_fraction_perm[x,])),
                  perm_sd=sapply(1:101, function(x) sd(signif_fraction_perm[x,])))

text_format <- element_text(face = "bold", color="black",size = 16)
cols <- c("Original"="black","Permuted"="steelblue")
p2=ggplot(data=ppa_df, aes(x=PPA, y=original)) + scale_x_reverse() + 
  geom_line(aes(colour="Original"), size=1.5) + theme_classic() + ylim(c(0,0.6)) +
  ylab("% regulatory variants (q<0.05)") + xlab("gPPA") + 
  geom_ribbon(aes(ymin = perm_mean-perm_sd, ymax = perm_mean+perm_sd), fill = "steelblue", alpha=0.5) + 
  geom_line(aes(y=perm_mean, colour="Permuted"), size=1.5) +
  theme(axis.text=text_format, axis.title=text_format, legend.text=text_format, legend.title = text_format, 
        legend.position=c(0.8,0.8)) + scale_colour_manual(name="P-values",values=cols) 
## add ribbon with the random permutaed data:
pdf("significant_variants_by_PPA.1000perms.cred_list.pdf")
p2
dev.off()

# pdf("significant_variants_by_PPA.1000perms.ADA.pdf")
# p2=ggplot(data=ppa_df, aes(x=PPA, y=original)) + scale_x_reverse() + 
#   geom_line(aes(colour="Original"), size=1.5) + theme_classic() + ylim(c(0,0.6)) +
#   ylab("% CNN regulatory variants") + xlab("gPPA") + 
#   geom_ribbon(aes(ymin = perm_mean-perm_sd, ymax = perm_mean+perm_sd), fill = "cadetblue3") + 
#   geom_line(aes(y=perm_mean, colour="Permuted"), size=1.5) +
#   theme(axis.text=text_format, axis.title=text_format, legend.text=text_format, legend.title = text_format, 
#         legend.position=c(0.8,0.8)) + scale_colour_manual(name="P-values",values=cols) 
# p2
# dev.off()

#### the same for FGWAS PPA:
#### going down in PPA from 1.00 to 0.01 in steps of 0.01 - what % of variants are significant <0.05 at this or higher PPA
signif_fraction_FGWAS=sapply(seq(from=1,to=0.00,by=-0.01), function(x) 
  length(which(subset(cred_df, FGWAS_PPA>=x & !is.na(FGWAS_PPA), lowest_Q)[,1]<0.05))/nrow(subset(cred_df, FGWAS_PPA>=x & !is.na(FGWAS_PPA))))
plot(1:101,signif_fraction_FGWAS, ylim=c(0,0.6), type="l")

AUC_fPPA=AUC(1:101, signif_fraction_FGWAS)

signif_fraction_perm_FGWAS=matrix(,nrow=101,ncol=n_perm)
for (i in 1:n_perm){
  signif_fraction_perm_FGWAS[,i]=sapply(seq(from=1,to=0.00,by=-0.01), function(x) 
    length(which(subset(cred_perm_df, FGWAS_PPA>=x,i+5)[,1]<0.05))/nrow(subset(cred_perm_df, FGWAS_PPA>=x)))
  print(i)
}

## p-value of this enrichment:
length(which(apply(signif_fraction_perm_FGWAS, 2, function(x) AUC(1:101,x)) >= AUC_fPPA))/n_perm


ppa_FGWAS_df=data.frame(PPA=seq(from=1,to=0.00,by=-0.01), original=signif_fraction_FGWAS,
                  perm_mean=sapply(1:101, function(x) mean(signif_fraction_perm_FGWAS[x,])),
                  perm_sd=sapply(1:101, function(x) sd(signif_fraction_perm_FGWAS[x,])))

text_format <- element_text(face = "bold", color="black",size = 16)
cols <- c("Original"="black","Permuted"="steelblue")
p6=ggplot(data=ppa_FGWAS_df, aes(x=PPA, y=original)) + scale_x_reverse() + 
  geom_line(aes(colour="Original"),size=1.5) + theme_classic() + ylim(c(0,0.6)) +
  ylab("% regulatory variants (q<0.05)") + xlab("fPPA") + 
  geom_ribbon(aes(ymin = perm_mean-perm_sd, ymax = perm_mean+perm_sd), fill = "steelblue", alpha=0.5) + 
  geom_line(aes(y=perm_mean, colour="Permuted"), size=1.5) +
  theme(axis.text=text_format, axis.title=text_format, legend.text=text_format, legend.title = text_format, 
        legend.position=c(0.8,0.8)) + scale_colour_manual(name="P-values",values=cols) 
p6
## add ribbon with the random permutaed data:
pdf("significant_variants_by_FGWAS_PPA.1000perms.list.pdf")
p6
dev.off()



#### going down in rank from 1.00 to 50 - what % of variants are significant <0.05 at this or higher PPA
### redone using the credible set structure, rather than subsetting data
signif_fraction_rank=rep(0,25)
for (i in 1:25){
  sig=length(which(unlist(lapply(cred, function(x) x[i,]$lowest_Q))<0.05))
  df=do.call(rbind,  lapply(cred, function(x) head(x, n=i)))
  sig=nrow(subset(df, lowest_Q<0.05))
  all=nrow(df)
  print(paste(i,sig,all))
  signif_fraction_rank[i]=sig/all
}
#plot(1:20,signif_fraction_rank[1:20], ylim=c(0.1,0.2), type="l")

signif_fraction_rank_perm=matrix(,nrow=25,ncol=n_perm)
for (i in 1:25){
  df=do.call(rbind,  lapply(cred_perm, function(x) head(x, n=i)))
  all=nrow(df)
  for (p in 1:n_perm){
    sig=length(which(df[,p+8]<0.05))
    #print(paste(p,i,sig,all))
    signif_fraction_rank_perm[i,p]=sig/all
  }
  print(i)
}
#plot(1:20,signif_fraction_rank[1:20], ylim=c(0.1,0.2), type="l")



rank_df=data.frame(rank=1:25, original=signif_fraction_rank,
                  perm_mean=sapply(1:25, function(x) mean(signif_fraction_rank_perm[x,])),
                  perm_sd=sapply(1:25, function(x) sd(signif_fraction_rank_perm[x,])))

text_format <- element_text(face = "bold", color="black",size = 16)
cols <- c("Original"="black","Permuted"="steelblue")
p3=ggplot(data=rank_df, aes(x=rank, y=original)) +  
  geom_line(aes(colour="Original")) + theme_classic() + ylim(c(0,0.2)) + xlim(c(1,25)) +
  ylab("% CNN regulatory variants") + xlab("Rank within credible sets") + 
  geom_ribbon(aes(ymin = perm_mean-perm_sd, ymax = perm_mean+perm_sd), fill = "steelblue", alpha=0.5) + 
  geom_line(aes(y=perm_mean, colour="Permuted")) +
  theme(axis.text=text_format, axis.title=text_format, legend.text=text_format, legend.title = text_format, 
        legend.position=c(0.8,0.2)) + scale_colour_manual(name="P-values",values=cols) 
p3
## add ribbon with the random permutaed data:
pdf("significant_variants_by_rank.1000perms.pdf")
p3
dev.off()

pdf("significant_variants_by_rank.1000perms.ADA.pdf")
p3=ggplot(data=rank_df, aes(x=rank, y=original)) +  
  geom_line(aes(colour="Original"), size=1.5) + theme_classic() + ylim(c(0,0.22)) + xlim(c(1,25)) +
  ylab("% CNN regulatory variants") + xlab("Rank within credible sets") + 
  geom_ribbon(aes(ymin = perm_mean-perm_sd, ymax = perm_mean+perm_sd), fill = "cadetblue3") + 
  geom_line(aes(y=perm_mean, colour="Permuted"), size=1.5) +
  theme(axis.text=text_format, axis.title=text_format, legend.text=text_format, legend.title = text_format, 
        legend.position=c(0.8,0.2)) + scale_colour_manual(name="P-values",values=cols) 
p3
dev.off()


### similar plots for just the IS vs IR loci

IS=scan("IS_loci.hard_cluster.txt", what="character")
IA=scan("IA_loci.hard_clustering.txt", what="character")

cred_IS=cred[subset(locus_ann, locus %in% IS, IndexSNP)[,1]]
cred_IA=cred[subset(locus_ann, locus %in% IA, IndexSNP)[,1]]

signif_fraction_IS=rep(0,25)
for (i in 1:25){
  df=do.call(rbind, lapply(cred_IS, function(x) head(x, n=i)))
  sig=nrow(subset(df, lowest_Q<0.05))
  all=nrow(df)
  print(paste(i, sig, all))
  signif_fraction_IS[i]=sig/all
}

signif_fraction_IA=rep(0,25)
for (i in 1:25){
  df=do.call(rbind, lapply(cred_IA, function(x) head(x, n=i)))
  sig=nrow(subset(df, lowest_Q<0.05))
  all=nrow(df)
  print(paste(i, sig, all))
  
  signif_fraction_IA[i]=sig/all
}

plot(1:20,signif_fraction_IS , type="l", ylim=c(0,0.25))
lines(1:20, signif_fraction_IA, col="grey")

# signif_fraction_rank_IS=sapply(seq(from=1,to=50,by=1), function(x) 
#   length(which(subset(data, locus %in% IS & credset_rank<=x, lowest_Q)[,1]<0.05))/nrow(subset(data, locus %in% IS & credset_rank<=x)))
# signif_fraction_rank_IA=sapply(seq(from=1,to=50,by=1), function(x) 
#   length(which(subset(data, locus %in% IA & credset_rank<=x, lowest_Q)[,1]<0.05))/nrow(subset(data, locus %in% IA & credset_rank<=x)))
# plot(1:50,signif_fraction_rank_IS , type="l", ylim=c(0,0.25), xlim=c(1,25))
# lines(1:50, signif_fraction_rank_IA, col="grey")

is_ia_df=data.frame(rank=1:50, IS=signif_fraction_IS, IA=signif_fraction_IA)
text_format <- element_text(face = "bold", color="black",size = 16)
cols <- c("Insulin Secretion (N=34)"="darkred","Insulin Action (N=23)"="steelblue")
p4=ggplot(data=is_ia_df, aes(x=rank, y=IS)) +  
  geom_line(aes(colour="Insulin Secretion (N=34)"), size=1.5) + theme_classic() + ylim(c(0,0.22)) + xlim(c(1,25)) +
  ylab("% CNN regulatory variants") + xlab("Rank within credible sets") + 
  geom_line(aes(y=IA, colour="Insulin Action (N=23)"), size=1.5) +
  theme(axis.text=text_format, axis.title=text_format, legend.text=text_format, legend.title = text_format, 
        legend.position=c(0.6,0.2)) + scale_colour_manual(name="Mechanism",values=cols) 
p4
## add ribbon with the random permutaed data:
pdf("IA_vs_IS_loci.rank.pdf")
p4
dev.off()


rank_df_merged=merge(rank_df, is_ia_df,by="rank")
cols <- c("Insulin Secretion (N=34)"="red","Insulin Action (N=23)"="purple","All signals (N=380)"="black","Permuted"="steelblue")
p5=ggplot(data=rank_df_merged, aes(x=rank, y=IS)) +  
  geom_ribbon(aes(ymin = perm_mean-perm_sd, ymax = perm_mean+perm_sd), fill = "steelblue", alpha=0.5) + 
  geom_line(aes(colour="Insulin Secretion (N=34)"), size=1.5) + theme_classic() + ylim(c(0,0.22)) + xlim(c(1,25)) +
  ylab("% regulatory variants (q<0.05)") + xlab("Rank within credible sets") + 
  geom_line(aes(y=perm_mean, colour="Permuted"), size=1.5) +
  geom_line(aes(y=original, colour="All signals (N=380)"), size=1.5) +
  geom_line(aes(y=IA, colour="Insulin Action (N=23)"), size=1.5) +
  theme(axis.text=text_format, axis.title=text_format, legend.text=text_format, legend.title = text_format, 
        legend.position=c(0.5,0.2)) + scale_colour_manual(name="",values=cols) 
p5
pdf("significant_variants_by_rank.IS_IA.1000perms.pdf")
p5
dev.off()



### compare also DeepSEA predictions - do we get better convergence using a specific tissue?
## read in all DeepSEA results - these are the multi-tissue predictions!!!
setwd("DeepSEA_HRC_results/")
deepsea_files=list.files(".",pattern="infile.vcf.out.funsig", recursive=T)

deepsea_all_res=data.frame()
for (I in deepsea_files){
  tmp=read.csv(I,h=T)
  deepsea_all_res=rbind(deepsea_all_res,tmp)
}
deepsea_all_res=deepsea_all_res[which(!duplicated(deepsea_all_res$name)),]
rownames(deepsea_all_res)=deepsea_all_res$name
names(deepsea_all_res)[4]="variant"
merged=merge(data, deepsea_all_res,by="variant")

### HERE #### append DeepSEA scores to cred list
for (sign in 1:length(cred)){
  cred[[sign]]$DeepSEA=rep(1,nrow(cred[[sign]]))

  for (var in 1:nrow(cred[[sign]])){
    cred[[sign]]$DeepSEA[var]=subset(deepsea_all_res, variant==cred[[sign]][var,]$variant, DeepSEA.score)[1,]

  }
  print(sign)
}
cred_df=rbindlist(cred)
signif_fraction_DeepSEA=sapply(seq(from=1,to=0.00,by=-0.01), function(x) 
  length(which(subset(cred_df, PPAg>=x, DeepSEA)[,1]<0.05))/nrow(subset(cred_df, PPAg>=x)))
### find random distribution:
### permute the DeepSEA scores 1000 times
### HERE ### 

perm_data_DeepSEA=subset(merged, , c(pos.x, lowest_Q, Genetic_PPA, FGWAS_PPA, DeepSEA.score))
n_perm=1000
set.seed(12345)
for (i in 1:n_perm){
    perm_data_DeepSEA[,5+i]=sample(perm_data_DeepSEA$DeepSEA.score)
    print(i)
}

cred_DeepSEA_perm=cred
for (c in 1:length(cred)){
  print(c)
  
  ### append the CNN permuted Q-values
  tmp=matrix(,nrow=nrow(cred[[c]]), ncol=n_perm)  
  for (i in 1:nrow(cred[[c]])){
    loc=paste(cred[[c]][i,]$Chr, cred[[c]][i,]$Pos,sep=":")
    tmp[i,]=as.numeric(subset(perm_data_DeepSEA, pos.x==loc )[1,c(6:(n_perm+5))])
  }
  cred_DeepSEA_perm[[c]]=cbind(cred[[c]], tmp)
}


## find what % of variants are regulatory at different PPA cut-offs
cred_DeepSEA_perm_df=rbindlist(cred_DeepSEA_perm)
cred_DeepSEA_perm_df=data.frame(cred_DeepSEA_perm_df)
#cred_perm_df=cred_perm_df[!(duplicated(paste0(cred_perm_df$Chr,":",cred_perm_df$Pos))),]

signif_fraction_DeepSEA_perm=matrix(,nrow=101,ncol=n_perm)
for (i in 1:n_perm){
  signif_fraction_DeepSEA_perm[,i]=sapply(seq(from=1,to=0.00,by=-0.01), function(x) 
    length(which(subset(cred_DeepSEA_perm_df, PPAg>=x,i+5)[,1]<0.05))/nrow(subset(cred_DeepSEA_perm_df, PPAg>=x)))
  print(i)
}

deepsea_ppa_df=data.frame(PPA=seq(from=1,to=0.00,by=-0.01), original=signif_fraction_DeepSEA,
                  perm_mean=sapply(1:101, function(x) mean(signif_fraction_DeepSEA_perm[x,])),
                  perm_sd=sapply(1:101, function(x) sd(signif_fraction_DeepSEA_perm[x,])))

text_format <- element_text(face = "bold", color="black",size = 16)
cols <- c("Original"="black","Permuted"="steelblue")
p10=ggplot(data=deepsea_ppa_df, aes(x=PPA, y=original)) + scale_x_reverse() + 
  geom_line(aes(colour="Original"), size=1.5) + theme_classic() + ylim(c(0,0.6)) +
  ylab("% CNN regulatory variants") + xlab("Genetic PPA") + 
  geom_ribbon(aes(ymin = perm_mean-perm_sd, ymax = perm_mean+perm_sd), fill = "steelblue", alpha=0.5) + 
  geom_line(aes(y=perm_mean, colour="Permuted"), size=1.5) +
  theme(axis.text=text_format, axis.title=text_format, legend.text=text_format, legend.title = text_format, 
        legend.position=c(0.8,0.8)) + scale_colour_manual(name="P-values",values=cols) 
## add ribbon with the random permutaed data:
pdf("deepsea.PPA_enrichment.pdf")
p10
dev.off()


##### repeat with ranks:
deepsea_signif_fraction_rank=rep(0,25)
for (i in 1:25){
  df=do.call(rbind,  lapply(cred, function(x) head(x, n=i)))
  sig=nrow(subset(df, DeepSEA<0.05))
  all=nrow(df)
  print(paste(i,sig,all))
  deepsea_signif_fraction_rank[i]=sig/all
}
plot(1:25,deepsea_signif_fraction_rank[1:25], ylim=c(0.1,0.3), type="l")

deepsea_signif_fraction_rank_perm=matrix(,nrow=20,ncol=n_perm)
for (i in 1:20){
  df=do.call(rbind,  lapply(cred_DeepSEA_perm, function(x) head(x, n=i)))
  all=nrow(df)
  for (p in 1:n_perm){
    sig=length(which(df[,p+9]<0.05))
    #print(paste(p,i,sig,all))
    deepsea_signif_fraction_rank_perm[i,p]=sig/all
  }
  print(i)
}
plot(1:20,deepsea_signif_fraction_rank_perm[1:20], ylim=c(0.1,0.2), type="l")


deepsea_rank_df=data.frame(rank=1:20, original=deepsea_signif_fraction_rank,
                   perm_mean=sapply(1:20, function(x) mean(deepsea_signif_fraction_rank_perm[x,])),
                   perm_sd=sapply(1:20, function(x) sd(deepsea_signif_fraction_rank_perm[x,])))

text_format <- element_text(face = "bold", color="black",size = 16)
cols <- c("Original"="black","Permuted"="steelblue")
p11=ggplot(data=deepsea_rank_df, aes(x=rank, y=original)) +  
  geom_line(aes(colour="Original")) + theme_classic() + ylim(c(0,0.3)) + xlim(c(1,20)) +
  ylab("% DeepSEA regulatory variants") + xlab("Rank within credible sets") + 
  geom_ribbon(aes(ymin = perm_mean-perm_sd, ymax = perm_mean+perm_sd), fill = "steelblue", alpha=0.5) + 
  geom_line(aes(y=perm_mean, colour="Permuted")) +
  theme(axis.text=text_format, axis.title=text_format, legend.text=text_format, legend.title = text_format, 
        legend.position=c(0.8,0.2)) + scale_colour_manual(name="P-values",values=cols) 
p11

##### compare for IA and IS:
cred_IS=cred[subset(locus_ann, locus %in% IS, IndexSNP)[,1]]
cred_IA=cred[subset(locus_ann, locus %in% IA, IndexSNP)[,1]]

deepsea_signif_fraction_IS=rep(0,20)
for (i in 1:20){
  df=do.call(rbind, lapply(cred_IS, function(x) head(x, n=i)))
  sig=nrow(subset(df, DeepSEA<0.05))
  all=nrow(df)
  print(paste(i, sig, all))
  deepsea_signif_fraction_IS[i]=sig/all
}

deepsea_signif_fraction_IA=rep(0,20)
for (i in 1:20){
  df=do.call(rbind, lapply(cred_IA, function(x) head(x, n=i)))
  sig=nrow(subset(df, DeepSEA<0.05))
  all=nrow(df)
  print(paste(i, sig, all))
  
  deepsea_signif_fraction_IA[i]=sig/all
}

plot(1:20,deepsea_signif_fraction_IS , type="l", ylim=c(0,0.3))
lines(1:20, deepsea_signif_fraction_IA, col="grey")

is_ia_df_deepsea=data.frame(rank=1:20, IS=deepsea_signif_fraction_IS, IA=deepsea_signif_fraction_IA)
text_format <- element_text(face = "bold", color="black",size = 16)
cols <- c("Insulin Secretion (N=34)"="darkred","Insulin Action (N=23)"="steelblue")
p12=ggplot(data=is_ia_df_deepsea, aes(x=rank, y=IS)) +  
  geom_line(aes(colour="Insulin Secretion (N=34)"), size=1.5) + theme_classic() + ylim(c(0,0.31)) + xlim(c(1,20)) +
  ylab("% CNN regulatory variants") + xlab("Rank within credible sets") + 
  geom_line(aes(y=IA, colour="Insulin Action (N=23)"), size=1.5) +
  theme(axis.text=text_format, axis.title=text_format, legend.text=text_format, legend.title = text_format, 
        legend.position=c(0.6,0.2)) + scale_colour_manual(name="Mechanism",values=cols) 
p12
## add ribbon with the random permutaed data:
pdf("IA_vs_IS_loci.rank.deepsea.pdf")
p12
dev.off()


deepsea_rank_df_merged=merge(deepsea_rank_df, is_ia_df_deepsea,by="rank")
cols <- c("Insulin Secretion (N=34)"="red","Insulin Action (N=23)"="purple","All signals (N=380)"="black","Permuted"="steelblue")
p13=ggplot(data=deepsea_rank_df_merged, aes(x=rank, y=IS)) +  
  geom_ribbon(aes(ymin = perm_mean-perm_sd, ymax = perm_mean+perm_sd), fill = "steelblue", alpha=0.5) + 
  geom_line(aes(colour="Insulin Secretion (N=34)"), size=1.5) + theme_classic() + ylim(c(0,0.32)) + xlim(c(1,20)) +
  ylab("% CNN regulatory variants") + xlab("Rank within credible sets") + 
  geom_line(aes(y=perm_mean, colour="Permuted"), size=1.5) +
  geom_line(aes(y=original, colour="All signals (N=380)"), size=1.5) +
  geom_line(aes(y=IA, colour="Insulin Action (N=23)"), size=1.5) +
  theme(axis.text=text_format, axis.title=text_format, legend.text=text_format, legend.title = text_format, 
        legend.position=c(0.5,0.2)) + scale_colour_manual(name="",values=cols) 
p13
pdf("significant_variants_by_rank.IS_IA.1000perms.deepSEA.pdf")
p13
dev.off()




######################   repeat this for PanIslets E-values  ##############

setwd("DeepSEA_HRC_results/")
panislets_files=list.files(".",pattern="*PanIslets", recursive=T)

panislets_all_res=data.frame()
for (I in panislets_files){
  tmp=read.csv(I,h=T)
  panislets_all_res=rbind(panislets_all_res,tmp)
}
panislets_all_res=panislets_all_res[which(!duplicated(panislets_all_res$name)),]
rownames(panislets_all_res)=panislets_all_res$name
names(panislets_all_res)[3]="variant"
merged=merge(merged, panislets_all_res,by="variant")

#### append PanIslets scores to cred list
for (sign in 1:length(cred)){
  cred[[sign]]$PanIslets=rep(1,nrow(cred[[sign]]))
  
  for (var in 1:nrow(cred[[sign]])){
    cred[[sign]]$PanIslets[var]=subset(panislets_all_res, variant==cred[[sign]][var,]$variant, PanIslets.DNase.None)[1,]
    
  }
  print(sign)
}
cred_df=rbindlist(cred)
signif_fraction_PanIslets=sapply(seq(from=1,to=0.00,by=-0.01), function(x) 
  length(which(subset(cred_df, PPAg>=x, PanIslets)[,1]<0.05))/nrow(subset(cred_df, PPAg>=x)))
### find random distribution:
### permute the PanIslets scores 1000 times

perm_data_PanIslets=subset(merged, , c(pos.x, lowest_Q, Genetic_PPA, FGWAS_PPA, PanIslets.DNase.None))
n_perm=1000
set.seed(123456)
for (i in 1:n_perm){
  perm_data_PanIslets[,5+i]=sample(perm_data_PanIslets$PanIslets.DNase.None)
  print(i)
}

cred_PanIslets_perm=cred
for (c in 1:length(cred)){
  print(c)
  
  ### append the CNN permuted Q-values
  tmp=matrix(,nrow=nrow(cred[[c]]), ncol=n_perm)  
  for (i in 1:nrow(cred[[c]])){
    loc=paste(cred[[c]][i,]$Chr, cred[[c]][i,]$Pos,sep=":")
    tmp[i,]=as.numeric(subset(perm_data_PanIslets, pos.x==loc )[1,c(6:(n_perm+5))])
  }
  cred_PanIslets_perm[[c]]=cbind(cred[[c]], tmp)
}


##### repeat with ranks:
panislets_signif_fraction_rank=rep(0,20)
for (i in 1:20){
  df=do.call(rbind,  lapply(cred, function(x) head(x, n=i)))
  sig=nrow(subset(df, PanIslets<0.05))
  all=nrow(df)
  print(paste(i,sig,all))
  panislets_signif_fraction_rank[i]=sig/all
}
plot(1:20,panislets_signif_fraction_rank[1:20], ylim=c(0,0.3), type="l")

panislets_signif_fraction_rank_perm=matrix(,nrow=20,ncol=n_perm)
for (i in 1:20){
  df=do.call(rbind,  lapply(cred_PanIslets_perm, function(x) head(x, n=i)))
  all=nrow(df)
  for (p in 1:n_perm){
    sig=length(which(df[,p+9]<0.05))
    #print(paste(p,i,sig,all))
    panislets_signif_fraction_rank_perm[i,p]=sig/all
  }
  print(i)
}
plot(1:20,panislets_signif_fraction_rank_perm[1:20], ylim=c(0,0.2), type="l")


panislets_rank_df=data.frame(rank=1:20, original=panislets_signif_fraction_rank,
                           perm_mean=sapply(1:20, function(x) mean(panislets_signif_fraction_rank_perm[x,])),
                           perm_sd=sapply(1:20, function(x) sd(panislets_signif_fraction_rank_perm[x,])))

text_format <- element_text(face = "bold", color="black",size = 16)
cols <- c("Original"="black","Permuted"="steelblue")
p21=ggplot(data=panislets_rank_df, aes(x=rank, y=original)) +  
  geom_line(aes(colour="Original")) + theme_classic() + ylim(c(0,0.15)) + xlim(c(1,20)) +
  ylab("% PanIslets (DeepSEA) regulatory variants") + xlab("Rank within credible sets") + 
  geom_ribbon(aes(ymin = perm_mean-perm_sd, ymax = perm_mean+perm_sd), fill = "steelblue", alpha=0.5) + 
  geom_line(aes(y=perm_mean, colour="Permuted")) +
  theme(axis.text=text_format, axis.title=text_format, legend.text=text_format, legend.title = text_format, 
        legend.position=c(0.8,0.2)) + scale_colour_manual(name="P-values",values=cols) 
p21

##### compare for IA and IS:
cred_IS=cred[subset(locus_ann, locus %in% IS, IndexSNP)[,1]]
cred_IA=cred[subset(locus_ann, locus %in% IA, IndexSNP)[,1]]

panislets_signif_fraction_IS=rep(0,20)
for (i in 1:20){
  df=do.call(rbind, lapply(cred_IS, function(x) head(x, n=i)))
  sig=nrow(subset(df, PanIslets<0.05))
  all=nrow(df)
  print(paste(i, sig, all))
  panislets_signif_fraction_IS[i]=sig/all
}

panislets_signif_fraction_IA=rep(0,20)
for (i in 1:20){
  df=do.call(rbind, lapply(cred_IA, function(x) head(x, n=i)))
  sig=nrow(subset(df, PanIslets<0.05))
  all=nrow(df)
  print(paste(i, sig, all))
  
  panislets_signif_fraction_IA[i]=sig/all
}

plot(1:20,panislets_signif_fraction_IS , type="l", ylim=c(0,0.3))
lines(1:20, panislets_signif_fraction_IA, col="grey")

is_ia_df_panislets=data.frame(rank=1:20, IS=panislets_signif_fraction_IS, IA=panislets_signif_fraction_IA)
text_format <- element_text(face = "bold", color="black",size = 16)
cols <- c("Insulin Secretion (N=34)"="darkred","Insulin Action (N=23)"="steelblue")
p22=ggplot(data=is_ia_df_panislets, aes(x=rank, y=IS)) +  
  geom_line(aes(colour="Insulin Secretion (N=34)"), size=1.5) + theme_classic() + ylim(c(0,0.2)) + xlim(c(1,20)) +
  ylab("% PanIslets regulatory variants") + xlab("Rank within credible sets") + 
  geom_line(aes(y=IA, colour="Insulin Action (N=23)"), size=1.5) +
  theme(axis.text=text_format, axis.title=text_format, legend.text=text_format, legend.title = text_format, 
        legend.position=c(0.6,0.2)) + scale_colour_manual(name="Mechanism",values=cols) 
p22
## add ribbon with the random permutaed data:
pdf("IA_vs_IS_loci.rank.PanIslets.pdf")
p22
dev.off()


panislets_rank_df_merged=merge(panislets_rank_df, is_ia_df_panislets,by="rank")
cols <- c("Insulin Secretion (N=34)"="red","Insulin Action (N=23)"="purple","All signals (N=380)"="black","Permuted"="steelblue")
p23=ggplot(data=panislets_rank_df_merged, aes(x=rank, y=IS)) +  
  geom_ribbon(aes(ymin = perm_mean-perm_sd, ymax = perm_mean+perm_sd), fill = "steelblue", alpha=0.5) + 
  geom_line(aes(colour="Insulin Secretion (N=34)"), size=1.5) + theme_classic() + ylim(c(0,0.2)) + xlim(c(1,20)) +
  ylab("% PanIslets regulatory variants") + xlab("Rank within credible sets") + 
  geom_line(aes(y=perm_mean, colour="Permuted"), size=1.5) +
  geom_line(aes(y=original, colour="All signals (N=380)"), size=1.5) +
  geom_line(aes(y=IA, colour="Insulin Action (N=23)"), size=1.5) +
  theme(axis.text=text_format, axis.title=text_format, legend.text=text_format, legend.title = text_format, 
        legend.position=c(0.5,0.2)) + scale_colour_manual(name="",values=cols) 
p23
pdf("significant_variants_by_rank.IS_IA.1000perms.PanIslets.pdf")
p23
dev.off()


## find what % of variants are regulatory at different PPA cut-offs
cred_PanIslets_perm_df=rbindlist(cred_PanIslets_perm)
cred_PanIslets_perm_df=data.frame(cred_PanIslets_perm_df)

signif_fraction_PanIslets_perm=matrix(,nrow=101,ncol=n_perm)
for (i in 1:n_perm){
  signif_fraction_PanIslets_perm[,i]=sapply(seq(from=1,to=0.00,by=-0.01), function(x) 
    length(which(subset(cred_PanIslets_perm_df, PPAg>=x,i+5)[,1]<0.05))/nrow(subset(cred_PanIslets_perm_df, PPAg>=x)))
  print(i)
}

panislets_ppa_df=data.frame(PPA=seq(from=1,to=0.00,by=-0.01), original=signif_fraction_PanIslets,
                            perm_mean=sapply(1:101, function(x) mean(cred_PanIslets_perm_df[x,])),
                            perm_sd=sapply(1:101, function(x) sd(cred_PanIslets_perm_df[x,])))

text_format <- element_text(face = "bold", color="black",size = 16)
cols <- c("Original"="black","Permuted"="steelblue")
p20=ggplot(data=panislets_ppa_df, aes(x=PPA, y=original)) + scale_x_reverse() + 
  geom_line(aes(colour="Original"), size=1.5) + theme_classic() + ylim(c(0,0.6)) +
  ylab("% CNN regulatory variants") + xlab("Genetic PPA") + 
  geom_ribbon(aes(ymin = perm_mean-perm_sd, ymax = perm_mean+perm_sd), fill = "steelblue", alpha=0.5) + 
  geom_line(aes(y=perm_mean, colour="Permuted"), size=1.5) +
  theme(axis.text=text_format, axis.title=text_format, legend.text=text_format, legend.title = text_format, 
        legend.position=c(0.8,0.8)) + scale_colour_manual(name="P-values",values=cols) 
## add ribbon with the random permutaed data:
pdf("PanIslets.PPA_enrichment.pdf")
p20
dev.off()












cred_perm_df$DeepSEA=cred_df$DeepSEA

signif_fraction_DeepSEA_perm=matrix(,nrow=101,ncol=n_perm)
for (i in 1:n_perm){
  signif_fraction_perm[,i]=sapply(seq(from=1,to=0.00,by=-0.01), function(x) 
    length(which(subset(cred_perm_df, PPAg>=x,i+5)[,1]<0.05))/nrow(subset(cred_perm_df, PPAg>=x)))
  print(i)
}
#### going down in PPA from 1.00 to 0.01 in steps of 0.01 - what % of variants are significant <0.05 at this or higher PPA
signif_fraction=sapply(seq(from=1,to=0.00,by=-0.01), function(x) 
  length(which(subset(data, Genetic_PPA>=x, lowest_Q)[,1]<0.05))/nrow(subset(data, Genetic_PPA>=x)))
plot(1:101,signif_fraction, ylim=c(0,0.6), type="l")


### make scatterplots of CNN q's and DeepSEA results


merged$deepsea_col = rep(rgb(0,0,0,50,maxColorValue=255),nrow(merged))
merged$deepsea_col[merged$lowest_Q<0.05 & merged$DeepSEA.score<0.05] = rgb(100,0,0,50,maxColorValue=255)
merged$deepsea_col[merged$lowest_Q<0.05 & merged$DeepSEA.score>=0.05] = rgb(0,100,0,50,maxColorValue=255)
merged$deepsea_col[merged$lowest_Q>=0.05 & merged$DeepSEA.score<0.05] = rgb(0,0,100,50,maxColorValue=255)

cor_p1 = ggplot(merged, aes(x=-log10(lowest_Q), y=-log10(DeepSEA.score))) + geom_point(color=merged$deepsea_col) +
  xlab("-log10(q-value) CNNs") +ylab("-log10(q-value) DeepSEA") + theme_bw() + 
   theme(axis.text=text_format, axis.title=text_format, legend.text=text_format, legend.title = text_format) +
  geom_text(x=250, y=5, label="r2=0.227")

merged$panislets_col = rep(rgb(0,0,0,50,maxColorValue=255),nrow(merged))
merged$panislets_col[merged$lowest_Q<0.05 & merged$PanIslets.DNase.None<0.05] = rgb(100,0,0,50,maxColorValue=255)
merged$panislets_col[merged$lowest_Q<0.05 & merged$PanIslets.DNase.None>=0.05] = rgb(0,100,0,50,maxColorValue=255)
merged$panislets_col[merged$lowest_Q>=0.05 & merged$PanIslets.DNase.None<0.05] = rgb(0,0,100,50,maxColorValue=255)

cor_p2 = ggplot(merged, aes(x=-log10(lowest_Q), y=-log10(PanIslets.DNase.None))) + geom_point(color=merged$panislets_col) +
  xlab("-log10(q-value) CNNs") +ylab("-log10(q-value) DeepSEA PanIslets") + theme_bw() + 
  theme(axis.text=text_format, axis.title=text_format, legend.text=text_format, legend.title = text_format) +
  geom_text(x=250, y=4.5, label="r2=0.223")




save.image("CNN_GWAS_convergence.Rdata")


###### compile figure - Figure 3 #######
setwd("~/Desktop/islets CNN/GWAS_CNN_convergence/")
text_format <- element_text( color="black",size = 16)
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

pdf("Fig3.pdf",width=12, height=5)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3))) # 3 rows, 5 columns
print(p2, vp = vplayout(1, 1))  # the big plot covers rows 1:2 and cols 1:3
print(p6, vp = vplayout(1, 2))
print(p5, vp = vplayout(1, 3))
dev.off()


###### compile Sup. FIgure 2 - comparison to DeepSEA
png("Sup.Fig2.png",width=12, height=12,units='in', res=300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2))) # 3 rows, 5 columns
print(cor_p1, vp = vplayout(1, 1))  
print(cor_p2, vp = vplayout(1, 2))
print(p13, vp = vplayout(2, 1))
print(p23, vp = vplayout(2, 2))
dev.off()

