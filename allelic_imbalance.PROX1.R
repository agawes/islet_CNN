prox_ase=read.table("../../luciferase_results/PROX1_ATAC_allele_count.txt",h=T,row.names = 1)
pdf("PROX_allelic_imbalance.pdf")
barplot(t(data.matrix(prox_ase[,c(2,1)])), main="",
        ylab="ATAC-seq reads", col=c("steelblue","red"), cex.axis=1.5, cex.lab=1.5,cex.names=1.5,
        legend = c("Major allele","Minor allele"), args.legend=list(bty="n", cex=1.5))
dev.off()

pdf("PROX_allelic_imbalance.h.pdf", height=3, width=6)
par(las=2, mar=c(5.1,8,1,1))
barplot(t(data.matrix(prox_ase[c(2,1),c(2,1)])), main="", horiz=TRUE, width=c(0.75, 0.75), ylim=c(0,2),
        xlab="ATAC-seq reads", col=c("grey","indianred"), cex.axis=1.5, cex.lab=1.5,cex.names=1.5,
        legend = c("Major allele","Minor allele"), args.legend=list(x="bottomright",bty="n", cex=1.5))
dev.off()



### allelic imbalance p-value:
rs177=binom.test(prox_ase["rs17712208","al1"], sum(prox_ase["rs17712208",]), (1/2), alternative="two.sided")$p.value
## p=2.92e-12
rs796=binom.test(prox_ase["rs79687284","al1"], sum(prox_ase["rs79687284",]), (1/2), alternative="two.sided")
## p=1.45e-04