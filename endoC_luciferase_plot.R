library(ggplot2)
library(dplyr)
library(ggsignif)

setwd("~/Desktop/islets CNN/luciferase_results/")
res=read.table("endoC_luciferase.txt",h=T)


summary <- res %>% # the names of the new data frame and the data frame to be summarised
  group_by(experiment) %>%   # the grouping variable
  summarise(mean = mean(luciferase),  # calculates the mean of each group
            sd = sd(luciferase), # calculates the standard deviation of each group
            n = n(),  # calculates the sample size per group
            SE = sd(luciferase)/sqrt(n())) # calculates the standard error of each group
 
summary$experiment=factor(summary$experiment, levels=c("GFP","EV","MTNR1B","WT","both_SNPs","rs17712208","rs79687284"))
levels(summary$experiment)[5]="both SNPs"
# Error bars represent standard error of the mean
text_format <- element_text(face = "bold", color="black",size = 16)

 ggplot(summary, aes(x=experiment, y=mean)) + 
   geom_bar(position=position_dodge(), stat="identity", fill=c("black","grey","white","grey","indianred","darkred","grey"),
            color=c("black"),width=0.75, size=1.2) +
   geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=0.2, size=1.2) + 
   ylab("Luciferase") + xlab("") + theme_classic() + 
   theme(axis.text.x=element_text(angle = 45, hjust = 1, face="bold",size=16),
         axis.text.y=text_format, axis.title=text_format)

 pdf("endoC_luciferase.for_paper.pdf")
 ggplot(summary, aes(x=experiment, y=mean)) + 
   geom_bar(position=position_dodge(), stat="identity", fill=c("white","white","grey","grey","indianred","darkred","grey"),
            color=c("black"),width=0.75, size=1.2) +
   geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=0.2, size=1.2) + 
   ylab("Luciferase") + xlab("") + theme_classic() + 
   theme(axis.text.x=element_text(angle = 45, hjust = 1, color="black",face="bold",size=16),
         axis.text.y=text_format, axis.title=text_format) +
   geom_signif(annotation="***", y_position=1.6, xmin=4, xmax=5,  tip_length = c(0.05, 0.4), textsize=6) +
   geom_signif(annotation="***", y_position=1.75, xmin=4, xmax=6,  tip_length = c(0.05, 0.7), textsize=6) +
   geom_signif(annotation="NS", y_position=1.9, xmin=4, xmax=7,  tip_length = c(0.05, 0.3), textsize=6) 
 dev.off()
 
 
 

 #### subset plot for presentation - skip MTNR1B
 pdf("endoC_luciferase.skip_controls.pdf", width=8,height=6.5)
 summary_subset=summary[c(1:3,5:7),]
 levels(summary_subset$experiment)[6]="rs17712208-A"
 levels(summary_subset$experiment)[7]="rs79687284-C"
 
 ggplot(summary_subset, aes(x=experiment, y=mean)) + 
   geom_bar(position=position_dodge(), stat="identity", fill=c("white","white","grey","indianred","indianred","grey"),
            color=c("black"),width=0.75, size=1.2) +
   geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=0.2, size=1.2) + 
   ylab("Relative Luciferase Units") + xlab("") + theme_classic() + 
   theme(axis.text.x=element_text(angle = 45, hjust = 1, color="black",face="bold",size=16),
         axis.text.y=text_format, axis.title=text_format) +
   geom_signif(annotation="***", y_position=1.6, xmin=3, xmax=4,  tip_length = c(0.05, 0.4), textsize=6) +
   geom_signif(annotation="***", y_position=1.75, xmin=3, xmax=5,  tip_length = c(0.05, 0.8), textsize=6) +
   geom_signif(annotation="NS", y_position=1.9, xmin=3, xmax=6,  tip_length = c(0.05, 0.3), textsize=6) 
dev.off()

### for paper, not bold
text_format <- element_text( color="black",size = 16)
pdf("endoC_luciferase.manuscript.pdf", width=6,height=6.5)
summary_subset=summary[c(1:3,5:7),]
levels(summary_subset$experiment)[6]="rs17712208-A"
levels(summary_subset$experiment)[7]="rs79687284-C"

ggplot(summary_subset, aes(x=experiment, y=mean)) + 
  geom_bar(position=position_dodge(), stat="identity", fill=c("white","white","grey","indianred","indianred","grey"),
           color=c("black"),width=0.7, size=1.2) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=0.2, size=1.2) + 
  ylab("Relative Luciferase Units") + xlab("") + theme_classic() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1, color="black",size=16),
        axis.text.y=text_format, axis.title=text_format) +
  geom_signif(annotation="***", y_position=1.6, xmin=3, xmax=4,  tip_length = c(0.05, 0.4), textsize=6) +
  geom_signif(annotation="***", y_position=1.75, xmin=3, xmax=5,  tip_length = c(0.05, 0.8), textsize=6) +
  geom_signif(annotation="NS", y_position=1.9, xmin=3, xmax=6,  tip_length = c(0.05, 0.3), textsize=6) 
dev.off()

