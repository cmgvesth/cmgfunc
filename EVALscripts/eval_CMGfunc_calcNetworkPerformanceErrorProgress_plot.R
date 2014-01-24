options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

Comment <- function(`@Comments`) {invisible()}

#################################### COMMENT START ####################################
Comment( `
#############################################################################
# Plot error progress for single network
#############################################################################

#--------------------------------------------------------------
# Input source
#--------------------------------------------------------------
python pre_CMGfunc_trainNetworks.py -i proteindata_clan_101.tab.C.train -t proteindata_clan_101.tab.C.test

#--------------------------------------------------------------
# Input
#--------------------------------------------------------------
head proteindata_clan_101.tab.A.train.progress
# Starting network training
Score:  0.128681682934  loop:  1
Score:  0.0941000672864  loop:  2
Score:  0.0779559119652  loop:  3
Score:  0.0704813414974  loop:  4
Score:  0.0657982863733  loop:  5
Score:  0.0626849095759  loop:  6

#--------------------------------------------------------------
# Output
#--------------------------------------------------------------
Scatter plot of MSE per iteration: plot_eval_CMGfunc_calcNetworkPerformanceErrorProgress.pdf

#--------------------------------------------------------------
# USAGE
#--------------------------------------------------------------
Rscript eval_CMGfunc_calcNetworkPerformanceErrorProgress_plot.R proteindata_clan_101.tab.A.train.progress

for x in *net.performance
do
Rscript eval_CMGfunc_calcNetworkPerformanceErrorProgress_plot.R $x
mv plot_eval_CMGfunc_calcNetworkPerformanceErrorProgress.pdf $x.plot_eval_CMGfunc_calcNetworkPerformanceErrorProgress.pdf
done
`)#################################### COMMENT END ####################################

string <- paste("plot_eval_CMGfunc_calcNetworkPerformanceErrorProgress.pdf")
pdf(string, width = 10, height=10)

library(ggplot2)
d <- read.table(args[1])

ggplot(d, aes(y= V2, x = seq(1, length(d$V2)))) + geom_point() +
geom_point(aes(y=d$V2), color="#0072B2") + 
scale_y_continuous(name="Network Mean squared error (MSE)")  + scale_x_continuous(name="Training data") +
geom_hline(yintercept=0.01, colour="red")

