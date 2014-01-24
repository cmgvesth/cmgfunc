options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

Comment <- function(`@Comments`) {invisible()}

#################################### COMMENT START ####################################
Comment( `
#############################################################################
# Histograms of MCC, PCC and MSE for a set of networks
#############################################################################

#--------------------------------------------------------------
# Input source
#--------------------------------------------------------------
python pre_CMGfunc_trainNetworks.py -i proteindata_clan_101.tab.C.train -t proteindata_clan_101.tab.C.test
grep -h MSE MAX/*progress* | grep test >> networkperformance.txt.raw
python eval_CMGfunc_calcNetworkPerformance.py -p networkperformance.txt.raw -s long > networkperformance.txt.tab

#--------------------------------------------------------------
# Input
#--------------------------------------------------------------
name	pcc	mcc	mse	sensitivity	specificity	precision
proteindata_clan_101. 	0.8654 	0.8201 	0.0363 	0.85243902439 	0.955568547208 	0.932414406403 	0.0444314527916 	0.85243902439
proteindata_clan_103. 	0.8093 	0.5717 	0.0313 	0.423924449108 	0.990353697749 	0.918181818182 	0.0096463022508 	0.423924449108
proteindata_clan_103. 	0.8733 	0.5334 	0.0313 	0.424973767051 	0.994908896034 	0.984602917342 	0.0050911039657 	0.424973767051

#--------------------------------------------------------------
# Output
#--------------------------------------------------------------
PDF of histograms: plot_eval_CMGfunc_calcNetworkPerformanceMultiNetHist_plot.pdf

#--------------------------------------------------------------
# USAGE
#--------------------------------------------------------------
Rscript eval_CMGfunc_calcNetworkPerformanceMultiNetHist_plot.R networkperformance.txt.tab
`) #################################### COMMENT END ####################################

#--------------------------------------------------------------
# Multiplot:
#---------------------------------------------------------------

#------ multiplot ------	# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Colors: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

library(ggplot2)
r <- read.table(args[1], header=T)


p1 <- ggplot(r, aes(x=pcc)) + geom_histogram(binwidth=0.01, aes(fill = ..count..), colour = "white") + 
scale_fill_gradient("Count", low = "white", high = "#009E73") + scale_x_continuous(name="PCC score, binwidth = 0.01") + 
scale_y_continuous(name="Nr. of networks") + geom_vline(xintercept=0.8, colour="red") + ggtitle("Pearson Correlation Coefficient histogram")

p2 <- ggplot(r, aes(x=mcc)) + geom_histogram(binwidth=0.01, aes(fill = ..count..), colour = "white") + 
scale_fill_gradient("Count", low = "white", high = "#0072B2") + scale_x_continuous(name="MCC score, binwidth = 0.01") + 
scale_y_continuous(name="Nr. of networks") + geom_vline(xintercept=0.8, colour="red") + ggtitle("Matthews Correlation Coefficient histogram")

p3 <- ggplot(r, aes(x=mse)) + geom_histogram(binwidth=0.001, aes(fill = ..count..), colour = "white") + 
scale_fill_gradient("Count", low = "white", high = "#CC79A7") + scale_x_continuous(name="MSE score, binwidth = 0.001") + 
scale_y_continuous(name="Nr. of networks") + geom_vline(xintercept=0.01, colour="red") + ggtitle("Mean Squared Error histogram")

p4 <- ggplot(r, aes(x=(sensitivity*100))) + geom_histogram(binwidth=1, aes(fill = ..count..), colour = "white") + 
scale_fill_gradient("Count", low = "white", high = "#E69F00") + scale_x_continuous(limits=c(10,100), name="Sensitivity percentage, binwidth = 1") + 
scale_y_continuous(name="Nr. of networks") + ggtitle("Sensitivity Percentage (TP/(TP+FN))*100 histogram")

p5 <- ggplot(r, aes(x=(specificity*100))) + geom_histogram(binwidth=1, aes(fill = ..count..), colour = "white") + 
scale_fill_gradient("Count", low = "white", high = "#56B4E9") + scale_x_continuous(limits=c(10,100), name="Specificity percentage, binwidth = 1") + 
scale_y_continuous(name="Nr. of networks") + ggtitle("Specificity Percentage (TN/(TN+PF))*100 histogram")

p6 <- ggplot(r, aes(x=(precision*100))) + geom_histogram(binwidth=1, aes(fill = ..count..), colour = "white") + 
scale_fill_gradient("Count", low = "white", high = "#D55E00") + scale_x_continuous(limits=c(10,100), name="precision percentage, binwidth = 1") + 
scale_y_continuous(name="Nr. of networks") + ggtitle("Precision Percentage (TP/(TP+FP))*100 histogram")

string <- paste("plot_eval_CMGfunc_calcNetworkPerformanceMultiNetHist_plot.pdf")
pdf(string, width = 50, height=30)
multiplot(p1, p4, p2,p6,p3,p5, cols=3)

string <- paste("plot_eval_CMGfunc_calcNetworkPerformance_PCC_Hist_plot.pdf")
pdf(string, width = 10, height=10)
plot(p1)
string <- paste("plot_eval_CMGfunc_calcNetworkPerformance_PCC_Hist_plot.png")
png(string)
plot(p1)

string <- paste("plot_eval_CMGfunc_calcNetworkPerformance_MCC_Hist_plot.pdf")
pdf(string, width = 10, height=10)
plot(p2)
string <- paste("plot_eval_CMGfunc_calcNetworkPerformance_MCC_Hist_plot.png")
png(string)
plot(p2)

string <- paste("plot_eval_CMGfunc_calcNetworkPerformance_MSE_Hist_plot.pdf")
pdf(string, width = 10, height=10)
plot(p3)
string <- paste("plot_eval_CMGfunc_calcNetworkPerformance_MSE_Hist_plot.png")
png(string)
plot(p3)

string <- paste("plot_eval_CMGfunc_calcNetworkPerformance_SEN_Hist_plot.pdf")
pdf(string, width = 10, height=10)
plot(p4)
string <- paste("plot_eval_CMGfunc_calcNetworkPerformance_SEN_Hist_plot.png")
png(string)
plot(p4)

string <- paste("plot_eval_CMGfunc_calcNetworkPerformance_SPE_Hist_plot.pdf")
pdf(string, width = 10, height=10)
plot(p5)
string <- paste("plot_eval_CMGfunc_calcNetworkPerformance_SPE_Hist_plot.png")
png(string)
plot(p5)

string <- paste("plot_eval_CMGfunc_calcNetworkPerformance_PRE_Hist_plot.pdf")
pdf(string, width = 10, height=10)
plot(p6)
string <- paste("plot_eval_CMGfunc_calcNetworkPerformance_PRE_Hist_plot.png")
png(string)
plot(p6)
