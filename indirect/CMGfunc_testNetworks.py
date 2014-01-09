#!/usr/bin/python

#------------------------ IMPORTS --------------------------------
from pybrain.datasets               import SupervisedDataSet #,ClassificationDataSet
from pybrain.datasets.classification    import SequenceClassificationDataSet
from pybrain.tools.validation       import testOnSequenceData
from pybrain.tools.customxml        import NetworkWriter, NetworkReader
from pybrain.tools.neuralnets       import NNregression, Trainer
from pybrain.tools.shortcuts        import buildNetwork
from pybrain.supervised.trainers    import BackpropTrainer #,RPropMinusTrainer
from pybrain.structure              import FeedForwardNetwork

from pylab          import ion, ioff, figure, draw, contourf, clf, show, hold, plot # pylab is needed for the graphical output
from scipy.stats    import pearsonr
from math           import sqrt
import os
from os import listdir
from os.path import isfile, join, isdir
import numpy
import argparse 
import logging 
import sys
import re
import subprocess

'''-----------------------------------------------------------------
#############################################################################
# Description
#############################################################################

#--------------------------------------------------------------
# Input source
#--------------------------------------------------------------
perl CMGfunc_normFeatures.pl -tab somefile.tab

OR FAKE:

for x in *A.test
do
LINE=$(grep -c "^1" $x | awk '{print $1*2}')
head -n $LINE $x > $x.pos
awk 'NR % 2 == 1' $x.pos > $x.tab
sed 's/^/name1\tname2\tname3\t/' $x.tab > $x.tab.vec
done

#--------------------------------------------------------------
# Input
#--------------------------------------------------------------
head trial.fsa.tab.vec
9272486_PF00005_PF00664_23	 trial.fsa	 9272486_PF00005_PF00664_23	 0	 0	 0	 1	 1	 0	 0	 1	 1	 0	 0	 0.232817869415808	 0.232817869415808	 0.00687285223367698	 0.0189003436426117	 0.0189003436426117	 0.00515463917525773	 0.122	 0.118	 0.143	 0.176	 0	 0	 0.243	 0.373	 0.134	 0.179	 0.148	 0.174	 0.024404	 0.373666666666667	 0.0166666666666667	 0.001	 0.001	 0.2571725	 0.084	 0.003	 0.052	 0.052	 0.040	 0.064	 0.015	 0.074	 0.043	 0.108	 0.046	 0.036	 0.021	 0.045	 0.062	 0.086	 0.058	 0.081	 0.009	 0.022	 0.000120274914089347	 0.0448594444444444	 0.0447405555555556	 0.903636363636364	 0.1	 0.170916	 0.553714285714286	 0.0136898461303436	 0.070	 0.103	 0.529	 0.120	 0.247	 0.333	 0.290	 0.206	 0.01455	 0.0321769759450172	 0.111218861852528	 0.130296425957675	 0.573875248096211	 0.175482871248489	 0.282611458472499	 
7242486_PF00005_PF08352_23	 trial.fsa	 7242486_PF00005_PF08352_23	 0.008	 0.211	 0.105	 0.788	 0.878	 0	 0.009	 0.788	 0.878	 0	 0	 0.475460122699386	 0.475460122699387	 0.00613496932515337	 0.0398773006134969	 0.0398773006134969	 0.00306748466257669	 0.152	 0.153	 0.114	 0.183	 0	 0	 0.175	 0.240	 0.102	 0.191	 0.125	 0.178	 1.2e-05	 0	 0	 0.001	 0.001	 0.2489275	 0.064	 0.028	 0.043	 0.071	 0.018	 0.080	 0.021	 0.074	 0.071	 0.110	 0.052	 0.018	 0.046	 0.064	 0.037	 0.058	 0.046	 0.074	 0.003	 0.021	 0.000131901840490798	 0.0281055555555556	 0.0271472222222222	 0.908181818181818	 0.1	 0.15794	 0.457428571428571	 0.00761666044285266	 0.043	 0.113	 0.549	 0.129	 0.209	 0.301	 0.298	 0.202	 0.00815	 0.0570398773006135	 0.068630698204604	 0.146175280858813	 0.586508469986988	 0.18215228566104	 0.248307075087616	

AND

Set of networkfiles
ll ../NETS/
total 5380
-rw-r--r-- 1 cmgfunc cmgfunc 45072 Nov 22 15:09 proteindata_clan_101.tab.A.net
-rw-r--r-- 1 cmgfunc cmgfunc 45006 Nov 22 15:09 proteindata_clan_103.tab.A.net
-rw-r--r-- 1 cmgfunc cmgfunc 45022 Nov 22 15:09 proteindata_clan_107.tab.A.net
-rw-r--r-- 1 cmgfunc cmgfunc 44992 Nov 22 15:09 proteindata_clan_109.tab.A.net
-rw-r--r-- 1 cmgfunc cmgfunc 44996 Nov 22 15:09 proteindata_clan_10.tab.A.net
....
#--------------------------------------------------------------
# Output
#--------------------------------------------------------------
head data_seq_trial_clan.fsa.tab.vec.res
Testseq:	10002486_PF02686_NULL 	Network:	proteindata_clan_157.tab.A.net 	Score:	0.87336622785
Testseq:	10002486_PF02686_NULL 	Network:	proteindata_clan_190.tab.A.net 	Score:	1.05753914864
Testseq:	10002486_PF02686_NULL 	Network:	proteindata_clan_214.tab.A.net 	Score:	0.85312832223
Testseq:	10002486_PF02686_NULL 	Network:	proteindata_clan_255.tab.A.net 	Score:	0.858665730197
Testseq:	10002486_PF02686_NULL 	Network:	proteindata_clan_509.tab.A.net 	Score:	0.983983750904


#--------------------------------------------------------------
# USAGE
#--------------------------------------------------------------
CMGfunc_testNetworks.py -i file.fsa.tab.vec -n /home/cmgfunc/CMGfunc/NETS/ -m 0.8

for x in NETS/proteindata_clan_* 
do 
echo $x
python pybrain_testNetworks.py -i data_seq_trial_clan.fsa.tab.vec -n $x
done
-----------------------------------------------------------------'''

#------------------------ OPTIONS --------------------------------
parser = argparse.ArgumentParser(description='# Get input and output file names.')
parser.add_argument('-i', required=True, help='Input filename <uarch10938.vec>')
parser.add_argument('-n', required=True, help='Path to network files')
parser.add_argument('-m', required=True, help='Score ')
args = parser.parse_args()
#-----------------------------------------------------------------

#.............. Output file ..............
out_file_name = args.i + ".res"
outall_file_name = args.i + ".res.all"

# If file already found, exit and ask to delete
if os.path.isfile(out_file_name):
    print "# WARNING CMGfunc_testNetworks.pl: resultfile already exists: ", out_file_name
    print "# WARNING CMGfunc_testNetworks.pl: will not re-calculate for sequences already in file, delete and re-run if desired"
    out_file = open(out_file_name, "a")
    outall_file = open(outall_file_name, "a")

else :
    out_file = open(out_file_name, "w")
    outall_file = open(outall_file_name, "w")

# Test that netpath exists
if not os.path.isdir(args.n):
    print "# ERROR: path to networks does not exist: ", args.n 
    sys.exit()

#.............. Network files ..............

# Get list of network files from a given directory
netfiles = [ f for f in listdir(args.n) if isfile(join(args.n,f)) ]	

#--------------------- READ TEST DATA ---------------------------
testvector = open( args.i )
line = testvector.readline()
seqfeats = []
i = 1



print >> out_file, "# Comparison started for file: ",args.i


while line:
    elements = line.split()
    seqname = elements[0]

    if os.path.isfile(out_file_name):
        out_file_read=open(out_file_name).read()
        if seqname in out_file_read:
           # print "Sequence already compared, %s, continue to next"%seqname
            line = testvector.readline()
            i=i+1
            continue

    seqfeats = map(float,elements[3:])

    max_score=0
    max_net = "none"
    for netfile in netfiles:

        net = NetworkReader.readFrom(args.n + netfile)
        score = net.activate(seqfeats)

        if (score > float(args.m)):
            print >> outall_file, "Testseq:\t", seqname, "\tNetwork:\t" , netfile ,"\tScore:\t", score        

        if (float(score) > max_score):
            max_score = float(score)
            max_net = netfile


    if (max_score > float(args.m)):
        print >> out_file, "Testseq:\t", seqname, "\tNetwork:\t" , max_net ,"\tScore:\t", max_score        
    else:
        print >> out_file, "Nomatch Testseq:\t", seqname, "\tNetwork:\t" , max_net ,"\tScore:\t", max_score        
       
    if (i % 100 == 0):
        print "# Running sequence number: ", i

    i=i+1
    line = testvector.readline()   

print >> out_file, "# Comparison ended for file: ",args.i
out_file.flush()

out_file.close()
outall_file.close()

