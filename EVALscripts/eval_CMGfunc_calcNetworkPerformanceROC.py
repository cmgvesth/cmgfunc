#!/usr/bin/python

#------------------------ IMPORTS --------------------------------
from pybrain.datasets             import ClassificationDataSet
from pybrain.utilities            import percentError
from pybrain.tools.shortcuts      import buildNetwork
from pybrain.supervised.trainers  import BackpropTrainer
from pybrain.structure.modules    import SoftmaxLayer
from pybrain.structure            import FeedForwardNetwork
from pybrain.tools.neuralnets     import NNregression, Trainer
from pylab                        import plot, hold, show
from scipy                        import sin, rand, arange
from pybrain.datasets.classification            import SequenceClassificationDataSet
from pybrain.structure.modules   import LSTMLayer, SoftmaxLayer
from pybrain.supervised          import RPropMinusTrainer
from pybrain.tools.validation    import testOnSequenceData
from pybrain.tools.shortcuts     import buildNetwork
from pybrain.datasets import SupervisedDataSet
from pybrain.tools.customxml     import NetworkWriter
from pybrain.tools.customxml     import NetworkReader
from pylab          import ion, ioff, figure, draw, contourf, clf, show, hold, plot # pylab is needed for the graphical output
from scipy          import diag, arange, meshgrid, where
from scipy.stats import pearsonr
import numpy
import argparse 
import logging 
from math import sqrt


'''
#############################################################################
# ROC curves
#############################################################################

#--------------------------------------------------------------
# Input source
#--------------------------------------------------------------
python pre_CMGfunc_trainNetworks.py -i proteindata_clan_101.tab.C.train -t proteindata_clan_101.tab.C.test -n networkperformance.txt

#--------------------------------------------------------------
# Input
#--------------------------------------------------------------
head h050_h150_100epo_rate01_setC/proteindata_clan_101.tab.C.net.performance
# Test data:proteindata_clan_101.tab.C.test
1.0 	0.989161342546
1.0 	0.994588789321
1.0 	0.935874858693
1.0 	0.988469570463
1.0 	0.991590932601
1.0 	1.02634602866
1.0 	1.08035848357
…….

#--------------------------------------------------------------
# Output
#--------------------------------------------------------------
#cut pcc mcc precision FPR TPR tp tn fp np rp rn
0.0 0.7705 0 1.0 0 0.961565668474 4078 0 0 163 4241 0
0.01 0.7705 0 0.211558307534 0.893306050862 1.0 820 365 3056 0 820 3421
0.02 0.7705 0 0.232756173716 0.790119847998 1.0 820 718 2703 0 820 3421
0.03 0.7705 0 0.262904777172 0.672025723473 1.0 820 1122 2299 0 820 3421
0.04 0.7705 0 0.292648108494 0.579362759427 1.0 820 1439 1982 0 820 3421
0.05 0.7705 0 0.323599052881 0.501023092663 1.0 820 1707 1714 0 820 3421
0.06 0.7705 0 0.348195329087 0.448699210757 1.0 820 1886 1535 0 820 3421
0.07 0.7705 0.4712 0.3733 0.4013 0.9976 818 2048 1373 2 820 3421
…..

#--------------------------------------------------------------
# USAGE
#--------------------------------------------------------------
python eval_CMGfunc_calcNetworkPerformanceROC.py -p file.net.performance

for x in *net.performance
do
python eval_CMGfunc_calcNetworkPerformanceROC.py -p $x
Rscript ../eval_CMGfunc_calcNetworkPerformanceROC_plot.R $x.ROC
mv plot_eval_calcNetworkPerformanceROC.pdf $x.plot_eval_calcNetworkPerformanceROC.pdf
done
'''

#-----------------------------------------------------------------

#------------------------ OPTIONS --------------------------------
parser = argparse.ArgumentParser(description='# Get input and output file names.')
parser.add_argument('-p', required=True, help='Network performance file: -p clan138.net.performance')
args = parser.parse_args()
#-----------------------------------------------------------------



#--------------------- READ TRAIN DATA ---------------------------
perfile = open( args.p )
perline = perfile.readline()
perline = perfile.readline()
predicted_num = []
target_num = []


while perline:

    elements = perline.split()
    predicted_num.append(float(elements[1]))
    target_num.append(float(elements[0]))
    perline = perfile.readline()

FHROC = open(args.p+".ROC", "w")

print >> FHROC,  "#cut pcc mcc precision FPR TPR tp tn fp np rp rn"
#for cut in [x * 0.01 for x in range(0, 10)]:
for cut in [x * 0.01 for x in range(0, 100)]:
    actual    = (numpy.array(target_num) >= float(cut))
    predicted = (numpy.array(predicted_num) >= float(cut))

    rp, rn, tp, fp , fn ,tn, pcc, mcc = 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0

    for i in range(len(actual)):
        if str(actual[i]) == "True": rp += 1
        if str(actual[i]) == "False": rn += 1
        if str(actual[i]) == "True" and str(predicted[i]) == "True": tp += 1 
        if str(actual[i]) == "True" and str(predicted[i]) == "False": fn += 1 
        if str(actual[i]) == "False" and str(predicted[i]) == "True": fp += 1 
        if str(actual[i]) == "False" and str(predicted[i]) == "False": tn += 1 
        #print tp, "\t", tn, "\t", fp, "\t", fn
        

    # Measures
    if (tp==0 or tn==0 or fp==0 or fn==0):
        mcc = 0
        if (tp==0 and fn==0): truePositiveRate = 0
        else: truePositiveRate = round(tp/(tp+fn), 4)

        if (tn==0 and fp==0): falsePositiveRate = 0
        else: falsePositiveRate = round(fp/(fp+tn),4)

        if (tp==0 and fp==0): precision = 0
        else: precision = round(tp/(tp+fp), 4)

    else:
        mcc = round( (tp*tn - fp*fn) / float(sqrt( (tp+fp) * (tp+fn) * (tn+fp) * (tn+fn) )), 4)
        precision = round(tp/(tp+fp), 4)
        falsePositiveRate = round(fp/(fp+tn), 4)
        truePositiveRate = round(tp/(tp+fn), 4)
        #print rp, rn, tp, fp , fn ,tn, round(mcc, 4), round(tp/fp, 4)
    
    pcc, PCCp_value = pearsonr(target_num,predicted_num) 
    
    print >> FHROC, cut,round(pcc, 4),mcc,precision, falsePositiveRate, truePositiveRate, int(tp), int(tn), int(fp), int(fn), int(rp), int(rn)
    #print cut,"\t",pcc,"\t",mcc,"\t",precision,"\t", falsePositiveRate, "\t", truePositiveRate, "\t",tp, "\t", tn, "\t", fp, "\t", fn, "\t", rp, "\t", rn 



