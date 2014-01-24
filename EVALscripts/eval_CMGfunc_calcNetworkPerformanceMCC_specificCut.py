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
# Calculate MCC at different cutoffs
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
.....

#--------------------------------------------------------------
# Output
#--------------------------------------------------------------
#	proteindata_clan_101.tab.C.net.performance 	PCC:	0.7705 	MCC_:	0.7643 	TP_:	731 	TN_:	3167 	FP_:	254 	FN_:	89 	RN:3421 	RP:	820
#--------------------------------------------------------------
# USAGE
#--------------------------------------------------------------
python eval_CMGfunc_calcNetworkPerformanceMCC_specificCut.py -p proteindata_clan_101.tab.C.net.performance -mcc 0.8

'''

#-----------------------------------------------------------------

#------------------------ OPTIONS --------------------------------
parser = argparse.ArgumentParser(description='# Get input and output file names.')
parser.add_argument('-p', required=True, help='Network performance file: -p clan138.net.performance')
parser.add_argument('-mcc', required=True, help='Cutoff for postitive netgative: -mcc 0.8')
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

actual    = (numpy.array(target_num) >= float(args.mcc))
predicted = (numpy.array(predicted_num) >= float(args.mcc))

rp, rn, tp, fp , fn ,tn = 0,0,0,0,0,0

for i in range(len(actual)):
    if str(actual[i]) == "True":
        rp += 1
    if str(actual[i]) == "False":
        rn += 1

    if str(actual[i]) == "True" and str(predicted[i]) == "True": 
        tp += 1 
    if str(actual[i]) =="True" and str(predicted[i]) =="False": 
        fn += 1 
    if str(actual[i]) =="False" and str(predicted[i]) =="True": 
        fp += 1 
    if str(actual[i]) =="False" and str(predicted[i]) =="False": 
        tn += 1 

# Measures
mcc = float(tp*tn - fp*fn) / float(sqrt( (tp+fp) * (tp+fn) * (tn+fp) * (tn+fn) ))
PCCr_row, PCCp_value = pearsonr(target_num,predicted_num) 

print "#\t%s" % args.p, "\tPCC:\t%.4f" % PCCr_row , "\tMCC_:\t%.4f" % mcc , "\tTP_:\t%.0f" % tp , "\tTN_:\t%.0f" % tn , "\tFP_:\t%.0f" % fp , "\tFN_:\t%.0f" % fn , "\tRN:\t%.0f" % rn , "\tRP:\t%.0f" % rp

