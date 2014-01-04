#! /usr/freeware/bin/gawk -f

# 
# This is a generic post-processor for raw HOW output. The input is pasted
# HOW output files. So far it has been shown to work with NetPhos 2.0,
# NetPhosK 1.0 and ProP 1.0b predictions.
# 
# External variables:
# 
# 	AA	amino acids to grep		e.g. S, ST, ...
# 	ID	prediction id			e.g. unsp, PKA, ...
# 	N	number of how.out files		e.g. 4, 5, ...
#	FL	flank to display (seq context)	e.g. 4
# 	FSC	final score calculation		e.g. ballot_1
# 						(if empty use average) 
# 
# VERSION:		 5 Oct 2004	launch, K. Rapacki
# 			 5 Apr 2005	'ballot_1' method added K. Rapacki
# 

BEGIN {

  if (FL!="") {				# prepare sequence flank
     nfl = split(FL,a,",");
     if (nfl==1) {
        lflank = substr("--------",1,FL);
	rflank = substr("--------",1,FL);
        lflank_len = FL;
	rflank_len = FL;
     }
     else {
        lflank = substr("--------",1,a[1]);
	rflank = substr("--------",1,a[2]);
        lflank_len = a[1];
	rflank_len = a[2];
     }   
  }

}

NF==0 {	act=0; }			# turn off printing

/SINGLE WINDOW OUTPUT ACTIVITIES:/ {	# new sequence

  getline <"in.tab";			# get seq name and seq

  seqname = $1;
  seq = (lflank $2 rflank);

  act++;				# turn on printing

}

$3 == "-" { next; }			# residue masked out

$2 ~ ("^[" AA "]$") {			# line matching amino acid

  # calculate final score from individual network scores
  scores=""; sum = 0;			# reset scores and sum

  if (FSC=="ballot_1") {			# find largest diff
     for (i=0; i<(N*6); i+=6) {
         scores = (scores "   " $(i+5)-$(i+6));
	 if (mabs($(i+5)-$(i+6))>sum) {
	    sum = mabs($(i+5)-$(i+6));
            sel = i;
	 }
     }
     sum = $(sel+5)-$(sel+6);
  }
  else {				# average (default)
     for (i=0; i<(N*6); i+=6) {
         scores = (scores "   " $(i+5));
         sum+=$(i+5);
     }
     sum/=N;
  }

  # print: seqname, position, residue
  printf("%-20s %5d %s   ",seqname,$1,$2);

  # print: sequence context
  if (FL!="")
     printf((substr(seq,$1,lflank_len+1) \
            a[3] \
	    substr(seq,$1+lflank_len+1,rflank_len)));

  # print: final score, prediction id and individual scores
  printf("   %5.3f   %-8s%s\n",sum,ID,scores);

}

function mabs(x) {

  if (x<0)
     return -x;
  else
     return x;

}
