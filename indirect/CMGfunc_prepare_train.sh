#!/bin/bash

:<<'COMMENT'
#############################################################################
# Create train and test files
#############################################################################

#--------------------------------------------------------------
# Input source
#--------------------------------------------------------------
perl CMGfunc_calcFeatures.pl -fasta genomeA.proteins.fsa

#--------------------------------------------------------------
# Input
#--------------------------------------------------------------
head proteindata_clan_101.tab 
10002112	PF00588_PF08032	101	0.25	0.926	0.25	0.024	0.25	0.001	0.25	0.926	0.25	0.001	0.048	0.996	0.996	0.004	0	0	0	0.103	0.109	0.117	0.128	00	0.152	0.204	0.123	0.136	0.111	0.122	0.000416	0	0	0.004	0.004	0.2331	0.104	0.012	0.06	0.064	0.016	0.076	0.032	0.052	0.044	0.092	0.036	0.064	0.036	0.02	0.064	0.052	0.048	0.092	0.012	0.024	0.000208	0.05695	0.0565333	0.907273	0.1	0.101632	0.4335	0.00581429	0.052	0.124	0.528	0.14	0.208	0.288	0.296	0.228	0.00625	0.073408	0.0861429	0.158762	0.553789	0.20569	0.245498
10004027	PF00588_PF08032	101	0.062	0.926	0.75	0.024	0.115	0.001	0.073	0.926	0.75	0.001	0.048	0.996	0.996	0.004	0	0	0	0.102	0.108	0.115	0.128	00	0.15	0.204	0.121	0.137	0.11	0.122	0.000402	0	0	0.004	0.004	0.2331	0.1	0.016	0.056	0.068	0.016	0.072	0.032	0.056	0.044	0.092	0.036	0.06	0.036	0.02	0.064	0.06	0.048	0.088	0.012	0.024	0.000208	0.0570889	0.0565333	0.907273	0.1	0.113872	0.433643	0.00582768	0.052	0.124	0.524	0.14	0.212	0.288	0.296	0.228	0.00625	0.07352	0.085945	0.158908	0.551461	0.205218	0.247776

#--------------------------------------------------------------
# Output
#--------------------------------------------------------------
*tab.A.test , *tab.A.train and prepare.log
Test file: 25% of positives and double number of negative examples
Train file: 75% of positives and double number of negative examples

==> description.local.test <==
==> description.local.train <==
0.006 0.004 0.178 0.982 0.816 0.001 0.001 0.982 0.816 0.001 0.012 0.455975 0.455975 0.0125786 0.0691824 0.0435109 0.00628931 0.129 0.108 0.094 0.086 0 0 0.104 0.116 0.082 0.069 0.105 0.097 0.004 0.004 0.269812 0.075 0.019 0.088 0.057 0.057 0.031 0.038 0.082 0.031 0.119 0.05 0.044 0.031 0.057 0.013 0.082 0.038 0.075 0.013 0 0.000433962 0.03909 0.0384344 0.903636 0.1 0.208852 0.309 0.00376565 0.069 0.145 0.553 0.082 0.22 0.346 0.302 0.189 0.003975 0.116805 0.106783 0.179604 0.624077 0.113271 0.243471
1
0.012 0.004 0.032 0.982 0.955 0.001 0.001 0.982 0.955 0.001 0.012 0.993711 0.993711 0.00628931 0 0 0 0.129 0.11 0.095 0.088 0 0 0.102 0.125 0.083 0.068 0.106 0.1 0.004 0.004 0.257547 0.075 0.019 0.082 0.057 0.057 0.044 0.025 0.063 0.031 0.126 0.05 0.038 0.038 0.069 0.025 0.075 0.038 0.075 0.013 0 0.000433962 0.03909 0.0384344 0.905455 0.1 0.216812 0.313071 0.0037578 0.069 0.138 0.56 0.082 0.22 0.333 0.308 0.195 0.003975 0.11673 0.107006 0.172464 0.625545 0.115658 0.24709
1
.......
0.012 0.004 0.032 0.982 0.955 0.001 0.001 0.982 0.955 0.001 0.012 0.993711 0.993711 0.00628931 0 0 0 0.129 0.108 0.094 0.086 0 0 0.104 0.116 0.082 0.069 0.105 0.097 0.004 0.004 0.281917 0.069 0.019 0.082 0.057 0.057 0.031 0.038 0.075 0.031 0.126 0.031 0.044 0.038 0.063 0.013 0.069 0.05 0.094 0.013 0 0.000433962 0.03909 0.0384344 0.902727 0.1 0.205556 0.312786 0.00375947 0.069 0.138 0.553 0.082 0.226 0.365 0.283 0.182 0.003975 0.117245 0.106959 0.172387 0.62114 0.113457 0.253703
0

#--------------------------------------------------------------
# USAGE
#--------------------------------------------------------------
pre_CMGfunc_prepare_train.sh description.fsa.tab.vec neg.tab.vec
COMMENT

pos=$1
neg=$2

rm $pos.train $pos.test

total_pos=$(grep -c "fsa" $pos | cut -f1)

if [ $total_pos -gt 50000 ]; then
   echo "# To many seqs in $pos, nr. lines $total_pos" >> prepare.log 
   exit;
fi
if [ $total_pos -lt 50 ]; then
   echo "# To few seqs in $pos, nr. lines $total_pos" >> prepare.log 
   exit;
fi

#---------------------------------------------------------------------
#----------------- Calculate number of samples used for training
#---------------------------------------------------------------------
train_pos=$(echo "$total_pos*0.75" | bc | awk '{print int($1)}')	# 75% of the pos data is used for training
train_neg=$(echo "$total_pos*1.5" | bc | awk '{print int($1)}')		# 150% of the size of the pos set is used as neg data in training

#---------------------------------------------------------------------
#----------------- Calculate number os samples used for testing
#---------------------------------------------------------------------
test_pos=$(echo "$total_pos*0.25" | bc | awk '{print int($1)}')	# 25% of the pos data is used for training	
test_neg=$(echo "$total_pos*0.5" | bc | awk '{print int($1)}')	# 50% of the size of the pos set is used as neg data in test

#echo $total_pos, $test_pos, $test_neg,  $train_pos, $train_neg

#---------------------------------------------------------------------
#----------------- Tranining data
#---------------------------------------------------------------------

#----------------- Positive data -----------------
#echo "# Running positive examples for $pos.train"

shuf -n $train_pos $pos | grep -v number | perl -lane 'print "@F[3..$#F]"; print "1"'  >> $pos.train

#echo "# Created positive examples for $pos.train"

#----------------- Negative data -----------------
#echo "# Running negative examples for $pos.train"

shuf -n $train_neg $neg | grep -v number | perl -lane 'print "@F[0..$#F]"; print "0"' >> $pos.train

#echo "# Created negative examples for $pos.train"

#---------------------------------------------------------------------
#----------------- Test data
#---------------------------------------------------------------------

#----------------- Positive data -----------------
#echo "# Running positive examples for $pos.test"

shuf -n $test_pos $pos | grep -v number | perl -lane 'print "@F[3..$#F]"; print "1"'  >> $pos.test
#shuf -n $test_pos $x | grep -v number | perl -lane 'print "@F[3..$#F] "; print "1"'  >> $pos.test

#echo "# Created positive examples for $pos.test"

#----------------- Negative data -----------------
#echo "# Running negative examples for $pos.test"

shuf -n $test_neg $neg | grep -v number | perl -lane 'print "@F[0..$#F]"; print "0"' >> $pos.test

#echo "# Created negative examples for $pos.test"

echo $pos, "tot pos. ", $total_pos, "test pos. ", $test_pos, "test neg. ", $test_neg,  "train pos. ", $train_pos, "train neg. ", $train_neg, "test lines neg. " >> prepare.log

#---------------------------------------------------------------------
#----------------- Pattern replacement
#---------------------------------------------------------------------

#echo "# Running pattern replacement for $pos.train"
sed -i '/^[ \t]*$/d' $pos.train

#echo "# Running pattern replacement for $pos.test"
sed -i '/^[ \t]*$/d' $pos.test

