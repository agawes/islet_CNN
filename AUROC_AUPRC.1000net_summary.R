act=read.table("data/learn_islets_act.txt")
act=act[grep("chr2:",rownames(act)),]

features=read.table("data/samples.txt")
features=as.character(features$V1)

library("PRROC")
library(calibrate)

files<-list.files(path=".", patter=".validation.txt$",recursive=TRUE)

cnn_res = matrix(,nrow=1000, ncol=60)
for(i in 1:length(files)){
    df = read.table(files[i])
    for(f in 1:length(features)){
		roc<-roc.curve(scores.class0 = df[,f], weights.class0 = act[,f])
		pr<-pr.curve(scores.class0 = df[,f], weights.class0 = act[,f])
    	cnn_res[i,2*f-1]= roc$auc
		cnn_res[i,2*f]= pr$auc.integral
	}
}
rownames(cnn_res)=gsub("runs/validation_set_predictions/iter","", gsub(".validation.txt", "",files))
cnn_res=data.frame(cnn_res)
colnames(cnn_res)[seq(from=1,to=59,by=2)]=paste0(features,"_AUROC")
colnames(cnn_res)[seq(from=2,to=60,by=2)]=paste0(features,"_AUPRC")

write.table(cnn_res, file="1000nets.AUC_ROC_PR.chr2.txt",sep="\t",quote=F)


# find mean, min, max - ROC and PR AUC
features=read.table("data/samples.txt")
features=as.character(features$V1)
roc_AUC=data.frame(features=features, mean_ROC_AUC=apply(cnn_res, 2, mean)[seq(from=1,to=59,by=2)], min_ROC_AUC=apply(cnn_res, 2, min)[seq(from=1,to=59,by=2)], max_ROC_AUC=apply(cnn_res, 2, max)[seq(from=1,to=59,by=2)], sd_ROC_AUC=apply(cnn_res, 2, sd)[seq(from=1,to=59,by=2)])
pr_AUC=data.frame(features=features, mean_ROC_AUC=apply(cnn_res, 2, mean)[seq(from=2,to=60,by=2)], min_ROC_AUC=apply(cnn_res, 2, min)[seq(from=2,to=60,by=2)], max_ROC_AUC=apply(cnn_res, 2, max)[seq(from=2,to=60,by=2)], sd_ROC_AUC=apply(cnn_res, 2, sd)[seq(from=2,to=60,by=2)])

write.table(cbind(roc_AUC, pr_AUC), file="1000nets.AUC_summary.txt", sep="\t",quote=F,row.names=F)


## find the best overall performing network and plot the representative ROC and PR curves for each feature
head(sort(apply(cnn_res,1,mean),decreasing=T))
## best in terms of PR:
head(sort(apply(cnn_res[,seq(2,60,2)],1,mean),decreasing=T))
#       578       619       583       250       686       570
# 0.4652003 0.4651296 0.4649779 0.4649697 0.4645486 0.4645304

files[578]
# [1] "fs21_runs/validation_set_predictions/iter79.validation.txt"

iter1=read.table("fs21_runs/validation_set_predictions/iter79.validation.txt")
pdf("fs21_iter79.roc_pr.pdf")
for(f in 1:length(features)){
	roc<-roc.curve(scores.class0 = iter1[,f], weights.class0 = act[,f], curve=T, max.compute = T, min.compute = T, rand.compute = T)
	pr<-pr.curve(scores.class0 = iter1[,f], weights.class0 = act[,f], curve=T, max.compute = T, min.compute = T, rand.compute = T)

	plot(roc, rand.plot = TRUE, auc.main=T, main=paste0(features[f]," ROC"))
	plot(pr, rand.plot = TRUE,auc.main=T, main=paste0(features[f]," PR"))

}
dev.off()
