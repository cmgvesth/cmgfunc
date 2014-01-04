#! /usr/freeware/bin/gawk -f

# /usr/cbs/bio/src/ape-1.0/disp/netphos.awk
# 
# This is NetPhos 3.1 display script for APE. It is also intended as a template
# script for other servers.
# 
# INPUT =======================================================================
# Output of 'nnhowplayer6' sorted by <seq_name>, <residue_num> and,
# inversely, <score_1> (in that order) e.g. lines like:
# 
# . S 5NTD_HUMAN 557 0.365573 0.632854    CKI
# 
# annotation, residue, seq_name, residue_num, score_1, score_2, label
# 
# OUTPUT ======================================================================
# Two output formats are currently supported:
# 
# 	ssa	Selective Sequence Annotation, default
# 	gff	General Feature Format
# 		(see http://www.sanger.ac.uk/Software/formats/GFF/)
# 
# Besides, graphics files in GIF, one per sequence, can be generated (see
# below).
# 
# ASSUMPTIONS =================================================================
# The files that have to exist etc.
# 
#	'in.tab'		exists in $cwd (APETMP)
# 	'<server-ver.>.gnu'	exists is $APE/graphics
# 
# COMMAND LINE VARIABLES ======================================================
# These variables can be passed in from the command line
# 
#    variable	meaning						default
# -----------------------------------------------------------------------------
# 	B	display only the best score for each residue	0
# 	C	cut-off for displaying scores			0
# 	E	economy of space: no separator lines		0
# 	F	output format					""
# 	G 	generate graphics				0
#	N	source-version					-
#	S	include sequence in the output			0
# -----------------------------------------------------------------------------
# 
# VERSION HISTORY
# 	2005 Aug 29	launch, K. Rapacki
# 

BEGIN {	# get sequences and lengths, prepare flanks
	flank="----";
	flen=length(flank);
	while ( "cat in.tab" | getline ) {
	      if (F=="gff")			# no context in GFF
	         seq[$1]=$2;
	      else				# residue context 
	         seq[$1]=(flank $2 flank);	
	      len[$1]=length($2);
	}

	# print overall GFF header
	if (F=="gff") {
	   print "##gff-version 2";
	   print "##source-version " N;
	   "date \"+%Y-%m-%d\"" | getline;
	   print "##date " $0;
	}
	else
	   S=1;		# SSA contains sequence by default

	# make sure data files for plotting exist
	if (G) system("touch {S,T,Y}.dat");

}

{
  if ($5<C) next;			# ignore scores lower than cut-off

  # first line (new sequence) #################################################
  if ($3!=seqname) {

     # generate graphics for previous sequence
     if (G && (seqname!="")) {
     	close("S.dat");close("T.dat");close("Y.dat");
     	cmd = ("sed 's/_LEN_/" len[seqname] "/;s/_NAME_/" seqname "/' ");
     	cmd = (cmd ENVIRON["APE"] "/graphics/" N ".gnu >" N ".gnu");
     	system(cmd);
     	cmd = (ENVIRON["GNUPLOT"] " " N ".gnu | " ENVIRON["PPM2GIF"]);
     	system(cmd ">" N "." seqname ".gif 2>/dev/null");
     }

     # print entry header, GFF
     if (F=="gff") {
        print "##Type Protein",$3;
	if (S)					# print current sequence
	   printseq();
	printf("# seqname            source        feature      ");
	print "start   end   score  N/A   ?";
	printf("# ----------------------------------------------");
	print "-----------------------------";
     }

     # print entry header, SSA
     else {
        if (seqname!="")			# print previous seq and ass
	   printseqass();
        for (i=1; i<=len[$3]; i++)		# reset assignment
            ass[i]=".";
        print ">"$3 "\t" len[$3],"amino acids";
	print "#\n#",N,"prediction results\n#";
        print "# Sequence\t\t   # x   Context     Score   Kinase    Answer";
	printf("# ----------------------------------------------");
	print "---------------------";
     }
     
     seqname=$3;
  }

  # any line ##################################################################
  if (B && ($4==prev)) next;		# ignore scores lower than the best

  if (F=="gff")
     printf("%-20s %s  phos-%-7s %5d %5d   %5.3f  . .  %s\n",
            seqname,N,$NF,$4,$4,$5,$5>=0.5?"YES":" . ");
  else {
     if (!E && prev && ($4!=prev)) print "#";
     context=substr(seq[seqname],$4,2*flen+1);
     printf("# %-20s %5d %s   %s   %5.3f   %-7s    %s\n",
            seqname,$4,$2,context,$5,$NF,$5>=0.5?"YES":" . ");
     if ($5>=0.5) ass[$4]=substr(seq[seqname],$4+flen,1);
  }
  prev=$4;

  # collect graphics data
  if (G) print $4,$5 >($2 ".dat");

}

###############################################################################
END {	if (seqname!="") {
	   if (F!="gff")			# print previous seq and ass
	      printseqass();
	   if (G) {				# generate graphics, last seq
              close("S.dat");close("T.dat");close("Y.dat");
              cmd = ("sed 's/_LEN_/" len[seqname] "/;s/_NAME_/" seqname "/' ");
              cmd = (cmd ENVIRON["APE"] "/graphics/" N ".gnu >" N ".gnu");
              system(cmd);
              cmd = (ENVIRON["GNUPLOT"] " " N ".gnu | " ENVIRON["PPM2GIF"]);
              system(cmd ">" N "." seqname ".gif 2>/dev/null");
           }

	}
}

# prints sequence and assignment, SSA style ###################################
function printseqass() {

  gsub("-","",seq[seqname]);
  print "#";
  for (i=1; i<=len[seqname]; i+=50)
      printf("    %-50s   # %6d\n",substr(seq[seqname],i,50),i+49);
  for (i=1; i<=len[seqname]; i++) {
      if (i % 50 == 1 ) printf("%1  ");
      printf(ass[i]);
      if (i % 50 == 0 ) printf("   # %6d\n",i);
  }
  if (i % 50 != 0 ) print "";

}

# prints sequence, GFF style ##################################################
function printseq() {

  print "##Protein",$3;
  for (i=1; i<=len[$3]; i+=60)
      printf("##%-60s\n",substr(seq[$3],i,60));
  print "##end-Protein";
}

# end of script ###############################################################
