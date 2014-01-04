#!/tools/opt/anaconda/bin/python2.7

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
import os
import re
import sys
'''

'''
#-----------------------------------
#----------------- SUB: testNetwork
#-----------------------------------
def testNetwork( testdata ):

    predicted_num = []
    target_num = []

    for i in xrange(len(testdata['input'])):
        seq = testdata['input'][i]
        target = testdata['target'][i]
        score = net2.activate(seq)

        predicted_num.append(float(score[0]))
        target_num.append(float(target))
           
    actual    = (numpy.array(target_num) >= 0.9)
    predicted = (numpy.array(predicted_num) >= 0.9)

    rp, rn, tp09, fp09 , fn09 ,tn09 = 0,0,0,0,0,0

    for i in range(len(actual)):
        if str(actual[i]) == "True":
            rp += 1
        if str(actual[i]) == "False":
            rn += 1

        if str(actual[i]) == "True" and str(predicted[i]) == "True": 
            tp09 += 1 
        if str(actual[i]) =="True" and str(predicted[i]) =="False": 
            fn09 += 1 
        if str(actual[i]) =="False" and str(predicted[i]) =="True": 
            fp09 += 1 
        if str(actual[i]) =="False" and str(predicted[i]) =="False": 
            tn09 += 1 
    
    return rp, rn, tp09, fp09 , fn09 ,tn09, target_num, predicted_num


#-----------------------------------------------------------------

#------------------------ OPTIONS --------------------------------
parser = argparse.ArgumentParser(description='# Get input and output file names.')
parser.add_argument('-i', required=True, help='Input filename <uarch10938.train>')
parser.add_argument('-t', required=True, help='Test filename <uarch10938.test>')
parser.add_argument('-c', required=False, help='Do not re-calc network', default='True')
parser.add_argument('-n', required=False, help='name for network performance file', default='networkperformance.txt')
args = parser.parse_args()
#-----------------------------------------------------------------

i_train = args.i
t_test = args.t
n_net = args.n



tmprootname = (args.i.split("/"))[-1]
rootname = tmprootname.replace('.train', '')
netname= args.i.replace('.train', '')
FH = open(netname+".progress", 'w')
print >> FH, "# Starting network training"
#--------------------- READ TRAIN DATA ---------------------------
trainvector = open( i_train )
trainline = trainvector.readline()
seqfeats, target_num , n = [], [], 0

traindata = SupervisedDataSet(75, 1) # create network
nr_pos = 0

while trainline:

    elements = trainline.split()
    if ((n % 2) == 0):
        seqfeats = map(float,elements)
        #print len(seqfeats)
    if ((n % 2) == 1):
        target_num = float(elements[0])
	traindata.addSample(seqfeats, target_num)
        nr_pos = nr_pos+1
    
    n += 1
    
    trainline = trainvector.readline()

#--------------------- READ TEST DATA ---------------------------
testvector = open( t_test )
testline = testvector.readline()
seqfeats, target_num , n = [], [], 0
testdata = SupervisedDataSet(75, 1) # create network

while testline:

    elements = testline.split()
    if ((n % 2) == 0):
        seqfeats = map(float,elements)
    if ((n % 2) == 1):
        target_num = float(elements[0])
        testdata.addSample(seqfeats, target_num)
	
    n += 1

    testline = testvector.readline()

#--------------------- PRINT INFO ---------------------------
print >> FH, "# Number of training patterns: ", len(traindata)
print >> FH, "# Number of testing patterns: ", len(testdata)
print >> FH, "# Input and output dimensions: ", traindata.indim, testdata.outdim

#--------------------- BUILD NETWORK ---------------------------
net = buildNetwork(traindata.indim, 30, 30, traindata.outdim)
net.randomize()

#--------------------- BUILD TRAINER ---------------------------
trainer = BackpropTrainer( net, dataset=traindata , learningrate=0.1)

#--------------------- TRAIN on trainingset ---------------------------
mse, score, i = 1,1,0

if (nr_pos < 1000):
    cut_score = 0.0001
else:
    cut_score = 0.001
        
if (args.c == "True"):
    while (i < 1000 and score > cut_score):
        score = trainer.train()
        i=i+1
        #print "Score: ", score, " loop: ", i
        print >> FH ,  "Score: ", score, " loop: ", i
        FH.flush()
    mse=score   # Store last MSE value as network MSE
    mse = round(mse, 4)
    netfilename = netname+".mse_"+str(mse)+".local.net"
    
    NetworkWriter.writeToFile(trainer.module, netfilename)   # Save network as XML file
    
#--------------------- TEST NETWORK on testset ---------------------------
net2 = NetworkReader.readFrom(netfilename)

rp, rn, tp09, fp09 , fn09 ,tn09 = 0,0,0,0,0,0
rp, rn, tp09, fp09 , fn09 ,tn09, target, predicted = testNetwork(testdata)

FHper = open(netfilename+".performance", "a")

print >> FHper ,  "# Test data:%s" % t_test

for x in xrange(len(target)):
    print >> FHper ,  target[x], "\t", predicted[x]

# Measures
mcc09 = float(tp09*tn09 - fp09*fn09) / float(sqrt( (tp09+fp09) * (tp09+fn09) * (tn09+fp09) * (tn09+fn09) ))
PCCr_row, PCCp_value = pearsonr(target,predicted) 

print >> FH , "#\t%s" % args.t, "\tPCC:\t%.4f" % PCCr_row , "\tMCC_09:\t%.4f" % mcc09 , "\tMSE:\t%.4f" % float(mse), "\tTP_09:\t%.0f" % tp09 , "\tTN_09:\t%.0f" % tn09 , "\tFP_09:\t%.0f" % fp09 , "\tFN_09:\t%.0f" % fn09 , "\tRN:\t%.0f" % rn , "\tRP:\t%.0f" % rp

#--------------------- TEST NETWORK on trainset ---------------------------
rp, rn, tp09, fp09 , fn09 ,tn09 = 0,0,0,0,0,0
rp, rn, tp09, fp09 , fn09 ,tn09, target, predicted = testNetwork(traindata)

# Measures
mcc09 = float(tp09*tn09 - fp09*fn09) / float(sqrt( (tp09+fp09) * (tp09+fn09) * (tn09+fp09) * (tn09+fn09) ))
PCCr_row, PCCp_value = pearsonr(target,predicted) 

print >> FH , "#\t%s" % args.i, "\tPCC:\t%.4f" % PCCr_row , "\tMCC_09:\t%.4f" % mcc09 , "\tMSE:\t%.4f" % float(mse), "\tTP_09:\t%.0f" % tp09 , "\tTN_09:\t%.0f" % tn09 , "\tFP_09:\t%.0f" % fp09 , "\tFN_09:\t%.0f" % fn09 , "\tRN:\t%.0f" % rn , "\tRP:\t%.0f" % rp

FH.flush()

# http://biomunky.wordpress.com/2010/03/17/hmmm-pybrains/
# http://stackoverflow.com/questions/8139822/how-to-load-training-data-in-pybrain/8143012#8143012
