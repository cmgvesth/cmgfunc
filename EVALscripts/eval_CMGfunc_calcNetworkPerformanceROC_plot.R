options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

library(ggplot2)

string <- paste("plot_eval_calcNetworkPerformanceROC.pdf")
pdf(string, width = 10, height=10)

s <- read.table(args[1])
ggplot(s, aes(x=V5, y=V6)) + geom_point() + 
scale_x_continuous(name="FPR, fp/(fp+tn)", limits=c(0,1)) + 
scale_y_continuous(name="TPR, tp/(tp+fn)", limits=c(0,1))



