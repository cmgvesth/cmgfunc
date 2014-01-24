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
# Analyze performance of different networks
# Takes the output of training networks
# Gives statistics for PCC, MCC and MSE for the set of networks
# Also calculates how many networks are above specific thresholds
# LONG mode gives statistics for individual networks
#############################################################################

#--------------------------------------------------------------
# Input source
#--------------------------------------------------------------
python pre_CMGfunc_trainNetworks.py -i proteindata_clan_101.tab.C.train -t proteindata_clan_101.tab.C.test
grep -h MSE MAX/*progress* | grep test >> networkperformance.txt.raw

#--------------------------------------------------------------
# Input
#--------------------------------------------------------------
head networkperformance.txt.raw
#	/home/people/tammi/projects/CMGfunc/SetClans/proteindata_clan_101.tab.C.test 	PCC:	0.8328 	MCC_09:	0.7854 	MSE:	0.0104 	TP_09:	588 	TN_09:	11326 	FP_09:	74 	FN_09:	232 	RN:	11400 	RP:	820
#	/home/people/tammi/projects/CMGfunc/SetClans/proteindata_clan_103.tab.B.test 	PCC:	0.8347 	MCC_09:	0.8004 	MSE:	0.0122 	TP_09:	753 	TN_09:	11258 	FP_09:	142 	FN_09:	200 	RN:	11400 	RP:	953
#	/home/people/tammi/projects/CMGfunc/SetClans/proteindata_clan_104.tab.C.test 	PCC:	0.8437 	MCC_09:	0.8032 	MSE:	0.0203 	TP_09:	2464 	TN_09:	12266 	FP_09:	274 	FN_09:	682 	RN:	12540 	RP:	3146
#	/home/people/tammi/projects/CMGfunc/SetClans/proteindata_clan_105.tab.C.test 	PCC:	0.8714 	MCC_09:	0.8204 	MSE:	0.0223 	TP_09:	7768 	TN_09:	39363 	FP_09:	537 	FN_09:	2231 	RN:	39900 	RP:	9999
#	/home/people/tammi/projects/CMGfunc/SetClans/proteindata_clan_106.tab.C.test 	PCC:	0.8862 	MCC_09:	0.8733 	MSE:	0.0161 	TP_09:	2347 	TN_09:	11183 	FP_09:	217 	FN_09:	324 	RN:	11400 	RP:	2671


#--------------------------------------------------------------
# Output
#--------------------------------------------------------------
python eval_CMGfunc_calcNetworkPerformance.py -p networkperformance.txt.raw -mcc 0.8 -pcc 0.8 -mse 0.01 -s short
# set	pcc_max	pcc_mean	pcc_min	mcc_max	mcc_mean	mcc_min	mse_max	mse_mean	mse_min	pcc_count	mcc_count	mse_count	sen_mean	spe_mean	pre_mean
networkperformance.txt.raw 	0.9674 	0.862417582418 	0.6317 	0.9629 	0.644781318681 	0.0263 	0.0677 	0.0300296703297 	0.0111 	37 	13 	0 	0.577682872484 	0.987428607073 	0.943610462839
....

python eval_CMGfunc_calcNetworkPerformance.py -p networkperformance.txt.raw -mcc 0.8 -pcc 0.8 -mse 0.01 -s long
# name	pcc	mcc	mse	sensitivity	specificity	precision
proteindata_clan_101. 	0.8654 	0.8201 	0.0363 	0.85243902439 	0.955568547208 	0.932414406403 	0.0444314527916 	0.85243902439
proteindata_clan_103. 	0.8093 	0.5717 	0.0313 	0.423924449108 	0.990353697749 	0.918181818182 	0.0096463022508 	0.423924449108
proteindata_clan_103. 	0.8733 	0.5334 	0.0313 	0.424973767051 	0.994908896034 	0.984602917342 	0.0050911039657 	0.424973767051
....

#--------------------------------------------------------------
# USAGE
#--------------------------------------------------------------
python eval_CMGfunc_calcNetworkPerformance.py -p networkperformance.txt.raw -s short
python eval_CMGfunc_calcNetworkPerformance.py -p networkperformance.txt.raw -s long


'''

#-----------------------------------------------------------------

#------------------------ OPTIONS --------------------------------
parser = argparse.ArgumentParser(description='# Get input and output file names.')
parser.add_argument('-p', required=True, help='Network performance file: -p networkperformance.txt.raw')
#parser.add_argument('-pcc', required=True, help='Cutoff for postitive negative: -pcc 0.8')
#parser.add_argument('-mcc', required=True, help='Cutoff for postitive negative: -mcc 0.8')
#parser.add_argument('-mse', required=True, help='Cutoff for postitive negative: -mse 0.03')
parser.add_argument('-s', required=True, help='per network (-s long) or summary (-s short)')
args = parser.parse_args()
#-----------------------------------------------------------------

#--------------------- READ TRAIN DATA ---------------------------
perfile = open( args.p )
perline = perfile.readline()
perline = perfile.readline()

pcc_list,mcc_list,mse_list = [],[],[]
tp_list,tn_list,fp_list,fn_list,rn_list,rp_list = [],[],[],[],[],[]
sen_list, spe_list, pre_list = [],[],[]

if (args.s == "long"):
    print "name\tpcc\tmcc\tmse\tsensitivity\tspecificity\tprecision\tfalsePositiveRate\ttruePositiveRate"

while perline:
    elements = perline.split()
    name = elements[1]
    pcc,mcc,mse = float(elements[3]), float(elements[5]), float(elements[7])
    tp,tn,fp,fn,rn,rp = float(elements[9]), float(elements[11]), float(elements[13]), float(elements[15]), float(elements[17]), float(elements[19])

    pcc_list.append(pcc)
    mcc_list.append(mcc)
    mse_list.append(mse)
    tp_list.append(tp)
    tn_list.append(tn)
    fp_list.append(fp)
    fn_list.append(fn)
    rn_list.append(rn)
    rp_list.append(rp)

    sensitivity = (tp/(tp+fn))
    specificity = (tn/(tn+fp))
    precision = (tp/(tp+fp))
    falsePositiveRate = fp/(fp+tn)
    truePositiveRate = tp/(tp+fn)

    sen_list.append(sensitivity)
    spe_list.append(specificity)
    pre_list.append(precision)

    if (args.s == "long"):
        print name,"\t",pcc,"\t",mcc,"\t",mse,"\t",sensitivity,"\t",specificity,"\t",precision,"\t", falsePositiveRate, "\t", truePositiveRate
    perline = perfile.readline()

if (args.s == "short"):

    pcc_count = sum(i > float(0.8) for i in pcc_list)
    mcc_count = sum(i > float(0.8) for i in mcc_list)
    mse_count = sum(i < float(0.03) for i in mse_list)

    print "set\tpcc_max\tpcc_mean\tpcc_min\tmcc_max\tmcc_mean\tmcc_min\tmse_max\tmse_mean\tmse_min\tpcc_count\tmcc_count\tmse_count\tsen_mean\tspe_mean\tpre_mean"
    print args.p,"\t",max(pcc_list),"\t",numpy.mean(pcc_list),"\t",min(pcc_list),"\t",max(mcc_list),"\t",numpy.mean(mcc_list),"\t", \
    min(mcc_list),"\t",max(mse_list),"\t",numpy.mean(mse_list),"\t",min(mse_list),"\t",pcc_count/2,"\t",mcc_count/2,"\t",mse_count/2,"\t", \
    numpy.mean(sen_list),"\t",numpy.mean(spe_list),"\t",numpy.mean(pre_list)



