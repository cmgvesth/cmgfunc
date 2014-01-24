options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

Comment <- function(`@Comments`) {invisible()}

#################################### COMMENT START ####################################
Comment( `
#############################################################################
# Plot target and predicted for single network
#############################################################################

#--------------------------------------------------------------
# Input source
#--------------------------------------------------------------
python pre_CMGfunc_trainNetworks.py -i proteindata_clan_101.tab.C.train -t proteindata_clan_101.tab.C.test

#--------------------------------------------------------------
# Input
#--------------------------------------------------------------
proteindata_clan_101.tab.C.train.progress
# Test data:proteindata_clan_101.tab.C.test
1.0 	0.989161342546
1.0 	0.994588789321
1.0 	0.935874858693
1.0 	0.988469570463
1.0 	0.991590932601
1.0 	1.02634602866
1.0 	1.08035848357
1.0 	1.05998085496

#--------------------------------------------------------------
# Output
#--------------------------------------------------------------
Scatter plot of target and predicted: plot_eval_CMGfunc_calcNetworkPerformanceTargetPredicted.pdf

#--------------------------------------------------------------
# USAGE
#--------------------------------------------------------------
Rscript eval_CMGfunc_calcNetworkPerformanceTargetPredicted_plot.R proteindata_clan_101.tab.C.net.progress

for x in *.progress
do
Rscript eval_CMGfunc_calcNetworkPerformanceTargetPredicted_plot.R $x
mv plot_eval_CMGfunc_calcNetworkPerformanceTargetPredicted.pdf $x.plot_eval_CMGfunc_calcNetworkPerformanceTargetPredicted.pdf
done
`)#################################### COMMENT END ####################################

string <- paste("plot_eval_CMGfunc_calcNetworkPerformanceTargetPredicted.pdf")
pdf(string, width = 10, height=10)

library(ggplot2)
d <- read.table(args[1])

ggplot(d, aes(y= V2, x = seq(1, length(d$V2)))) + geom_point() +
geom_point(aes(y=d$V2), color="#0072B2") + 
scale_y_continuous(name="Network Mean squared error (MSE)")  + scale_x_continuous(name="Training data") +
geom_hline(yintercept=0.01, colour="red")
