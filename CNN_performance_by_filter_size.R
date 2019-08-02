setwd("~/Desktop/islets CNN/1000nets_results/")

perf=read.table("1000nets.AUC_ROC_PR.chr2.txt",h=T,row.names=1)
perf$fs=as.numeric(gsub("fs","",gsub("_.+","",rownames(perf))))

pdf("CNN_performance_by_filter_size.pdf")
for (feature in colnames(perf)[1:60]){
  print(feature)
  fs_values=c(7,9,11,13,15,17,19,21,23,25)
  plot(fs_values,aggregate(perf[,feature], by=list(fs=perf$fs), FUN=mean)$x, type="b", ylab=feature)
  points(fs_values,aggregate(perf[,feature], by=list(fs=perf$fs), FUN=mean)$x,pch=16)
}
dev.off()

