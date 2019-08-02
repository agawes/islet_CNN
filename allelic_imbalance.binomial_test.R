aseq_res = list.files("./", recursive=T)
aseq=list()
for (f in aseq_res){
	base=gsub("_mapped.+","",(gsub("sorted.","",gsub("_mapped.PILEUP.ASEQ","",f))))
	print(base)
	aseq[[base]] = read.table(f, h=T,sep="\t")
}

#vcf=read.table("../../HRC.PPA_above0.5.ASEQ.vcf",sep="\t")
#vcf=read.table("../../HRC.PPA_0.1_0.5.ASEQ.vcf",sep="\t")
vcf=read.table("../../HRC.PPA_below_0.1.ASEQ.vcf",sep="\t")

snps = levels(vcf$V3)

snp_ASE = data.frame()

for (snp in snps){
	print(snp)
	aseq_tmp = lapply(aseq, function(x) x[which(x$dbsnp == snp & x$af >0.05 &  x$af<0.95),])
	aseq_tmp = aseq_tmp[lapply(aseq_tmp,nrow)>0]
	if(length(aseq_tmp)>0){
		df <- data.frame(matrix(unlist(aseq_tmp), nrow=length(aseq_tmp), byrow=T))
		names(df) = names(aseq_tmp[[1]])
		### binomial test per row, and total
		al1=names(sort(colSums(df[,7:10]))[4])
		al2=names(sort(colSums(df[,7:10]))[3])
		if (sum(df[,al2])>=5 & sum(df[,al1]) >=5) {
			p=binom.test(sum(df[,al2]), sum(df[,al1])+sum(df[,al2]), (1/2), alternative="two.sided")$p.value
			print(p)
			#snp_p[which(snps == snp)] = p
			snp_ASE=rbind(snp_ASE, data.frame(snp=snp, al1=sum(df[,al1]), al2=sum(df[,al2]), ratio=sum(df[,al2])/(sum(df[,al1])+sum(df[,al2])), p=p))
		} 
		if (sum(df[,al2])<5){
			print("Low cov")
			#snp_p[which(snps == snp)] = NA
			snp_ASE=rbind(snp_ASE, data.frame(snp=snp, al1=sum(df[,al1]), al2=sum(df[,al2]), ratio=sum(df[,al2])/(sum(df[,al1])+sum(df[,al2])), p=NA))
		}
	} 
	if(length(aseq_tmp)==0){print("No hets")
		#snp_p[which(snps == snp)] = NA
		snp_ASE=rbind(snp_ASE, data.frame(snp=snp, al1=NA, al2=NA, ratio=NA, p=NA))
	}
}
head(snp_ASE[order(snp_ASE$p),])
#write.table(snp_ASE, file="HRC.PPA_above0.5.allelic_imbalance.txt",sep="\t",quote=F)
#write.table(snp_ASE, file="HRC.PPA_0.1_0.5.allelic_imbalance.txt",sep="\t",quote=F)
write.table(snp_ASE, file="HRC.PPA_below_0.1.allelic_imbalance.txt",sep="\t",quote=F)

aaa=read.table("/well/got2d/agata/Basset/Full_islet_dataset/rank_product.1000nets.PPA_above_0.5.ann.txt",sep="\t")
aaa=read.table("/well/got2d/agata/Basset/Full_islet_dataset/rank_product.1000nets.PPA_0.1_0.5.ann.txt",sep="\t")

aaa$V1
rownames(snp_ASE) = snp_ASE$snp
snp_ASE[aaa$V1,]

