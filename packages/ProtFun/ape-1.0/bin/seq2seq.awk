#! /usr/freeware/bin/gawk -f


# Command line parameters:
# 	B		bottom: min seq length allowed
# 	OUTDIR		output directory
# 	OUTPUTFORMAT	output format
# 	MAXNAMELEN      max name length
# 	M		masking instructions: <aalist>,<range>
# 	N		flag: replace non-nucleotide symbols with this char
# 	P		flag: replace non-aa symbols with this char
# 	PRES		preserve do not alter names, include comments
# 	T		top: max seq length allowed
# 	TRUNC		truncation: if 0 ignore, if negative from the end
# 	V		verbosity level
# 	X		flag: do not tolerate non-standard symbols
# 
# Functions:
# 	consider()
# 	reset()
# 	report()
# 	printfasta(OUTDIR,len,name,seq)
# 	printhow(OUTDIR,len,name,seq,ass)
# 	printtab(OUTDIR,len,name,seq,ass)
# 	printact(OUTDIR,len,name,seq,ass)
# 

BEGIN {

  # allowed sequence symbols
  NON_NUCL_ALPH	= "[^ACGTU]";
  NON_PROT_ALPH	= "[^ACDEFGHIKLMNPQRSTVWY]";

  # prepare output directory
  if (OUTDIR!="/dev/stdout") {
     system("rm -rf " OUTDIR);
     system("mkdir " OUTDIR);
  }

  # masking
  maa = "_";			# default: no amino acid will match
  mbottom = 99999;		# default: no position above ;-)
  mtop = 0;			# default: no position below ;-)

  s1 = split(M,a,",");		# first split: aminoacids to mask out
  if (a[1]!="")
     maa = ("[" a[1] "]");

  if (s1>1) {			# second split: positions to mask out
     split(a[2],b,"-");
     mbottom = b[1];
     mtop = b[2];
     if (mbottom=="")
        mbottom=1;
     if (mtop=="")
        mtop=99999;
  }

  # truncation
  gsub("^[-][-]*","-",TRUNC);	# no extended regexp in AWK ...

}

# ban empty lines =============================================================
!NF { print "seq2seq: empty input line, not allowed" >"/dev/stderr"; exit; }

# HOW entry header ============================================================
/^[ 1-9]?[ 0-9][ 0-9][ 0-9][ 0-9][0-9] [^ \t]/ {

  consider();			# print previous entry if applicable

  reset();			# reset flags and entry data

  eformat="how";
  ec[eformat]++;		# increment input entry count

  len = $1;			# sequence length
  name = $2;			# sequence name
  
  comment = "";			# comment
  for (cw=3; cw<=NF; cw++)
      comment = (comment " " $cw);

  i = 0;			# running index

  ishowseq = 1;			# flag: picking up HOW sequence

  next;

}

# FASTA entry header ==========================================================
/^>/ {

  consider();			# print previous entry if applicable

  reset();			# reset flags and entry data

  eformat="fasta";
  ec[eformat]++;		# increment entry count

  len = 0;			# sequence length (not specified in FASTA)
  name = substr($1,2);		# sequence name

  if (name=="")			# special case: nameless FASTA entry
     name = ("seq." ec["how"]+ec["fasta"]);

  comment = "";			# comment
  for (cw=2; cw<=NF; cw++)
      comment = (comment " " $cw);

  i = 0;			# running index

  isfastaseq = 1;		# flag: picking up FASTA sequence

  next;

}

###############################################################################
corruption { next; }
###############################################################################

# picking up HOW sequence =====================================================
ishowseq {

  # get and check line length .................................................
  llen = length($1);		# sequence line length (not more than 80)
  if (llen>80) {
     if (V)
        print "seq2seq: HOW sequence line too long in " name "\"" \
		>"/dev/stderr";
     corruption++;
     next;
  }

  # convert to upper-case and replace non-standard symbols ....................
  $1 = toupper($1);
  if (N!="")
     nonstd += gsub(NON_NUCL_ALPH,N,$1);
  else if (P!="")
     nonstd += gsub(NON_PROT_ALPH,P,$1);

  # load sequence .............................................................
  for (j=1; j<=llen; j++)
      seq[++i] = substr($1,j,1);

  # check total sequence length so far ........................................
  if (i>len) {
     if (V)
        print "seq2seq: HOW sequence too long in \"" name "\"" \
		>"/dev/stderr";
     corruption++;
  }
  else if (i==len) {
     ishowseq = 0;
     ishowass = 1;
     seqlen = i;
     i = 0;
  }

  next;
}

# picking up HOW assignment ===================================================
ishowass {

  # get and check line length .................................................
  llen = length($1);		# assignment line length (not more than 80)
  if (llen>80) {
     if (V)
        print "seq2seq: HOW assignment line too long in \"" name "\"" \
		>"/dev/stderr";
     corruption++;
     next;
  }
  
  # load assignment ...........................................................
  for (j=1; j<=llen; j++)
      ass[++i] = substr($1,j,1);

  # check total assignment length so far ......................................
  if (i>len) {
     if (V)
        print "seq2seq: HOW assignment too long in \"" name "\"" \
		>"/dev/stderr";
     corruption++;
     next;
  }

  asslen = i;

  next;
}

# picking up FASTA sequence ===================================================
isfastaseq {

  # remove whitespace and get line length .....................................
  gsub("[ \t]","",$0); 
  llen = length($0);

  # convert to upper-case and replace non-standard symbols ....................
  $0 = toupper($0);
  if (N!="")
     nonstd += gsub(NON_NUCL_ALPH,N,$1);
  else if (P!="")
     nonstd += gsub(NON_PROT_ALPH,P,$1);

  # load sequence .............................................................
  for (j=1; j<=llen; j++)
      seq[++i] = substr($0,j,1);

  len=i;				# keep it updated ...
  seqlen = i;

  next;
}

# picking up TAB entry ========================================================
{
  consider();			# print previous entry if applicable

  reset();			# reset flags and entry data

  eformat="tab";
  ec[eformat]++;		# increment input entry count

  name = $1;			# sequence name
  seqlen = length($2);		# sequence length
  asslen = length($3);		# assignment length

  comment = "";			# comments not allowed in this format

  # syntax checks .............................................................
  if (!seqlen) {
     if (V)
        print "seq2seq: no sequence found in entry \"" name "\"" \
		>"/dev/stderr";
     corruption++;
     next;
  }

  # convert to upper-case and replace non-standard symbols ....................
  $2 = toupper($2);
  if (N!="")
     nonstd += gsub(NON_NUCL_ALPH,N,$2);
  else if (P!="")
     nonstd += gsub(NON_PROT_ALPH,P,$2);

  # pick up sequence and assignment ...........................................
  for (i=1; i<=seqlen; i++)
      seq[i]=substr($2,i,1);

  if (seqlen==asslen)			# assignment OK
     for (i=1; i<=seqlen; i++)
         ass[i]=substr($3,i,1);
  else					# no valid assignment, create empty
     for (i=1; i<=seqlen; i++)
         ass[i]=".";

  len = seqlen;
  asslen = len;

  next;
}

# =============================================================================
END {	consider(); 
        if (V>1)
	   report();
}


###############################################################################
function consider() {
###############################################################################

  if (eformat=="")			# non-existent previous entry
     return;

  if (corruption) {			# count corrupted entries
     cec[eformat]++;
     return 0;
  }
  else if (seqlen==0) {			# count zero length entries
     if (V)
        print "seq2seq: entry \"" name "\" has zero length" \
		>"/dev/stderr";
     zec[eformat]++;
     return 0;
  }
  else if (T && (seqlen>T)) {		# count too long entries
     if (V)
        print "seq2seq: entry \"" name "\" is too long" \
		>"/dev/stderr";
     Tec[eformat]++;
     return 0;
  }
  else if (seqlen<B) {			# count too short entries
     if (V)
        print "seq2seq: entry \"" name "\" is too short" \
		>"/dev/stderr";
     Bec[eformat]++;
     return 0;
  }

  # fix name ..................................................................
  origname = name;			# special characters
  if (PRES)
     namesub = 0;
  else
     namesub = gsub("[^A-Za-z0-9+,-._]","_",name);
  if (namesub) {
     if (V)
        #print "seq2seq: entry name altered from \"" origname "\" to \"" \
		#name "\"" >"/dev/stderr";
     mec[eformat]++;
  }

  if (MAXNAMELEN) {			# length
     if (length(name)>MAXNAMELEN) {
        name = substr(name,1,MAXNAMELEN);
        if (V)
           print "seq2seq: entry name truncated from \"" origname "\" to \"" \
		name "\"" >"/dev/stderr";
        tec[eformat]++;
     }
  }

  if (name in present) {		# non-unique name
     if (V)
        print "seq2seq: occurrence",present[name]+1,"of \"" name "\"," \
		" discarded"  >"/dev/stderr";
     dec[eformat]++;
     present[name]++;
     return 0;
  }
  present[name]=1;

  # end of name fixing ........................................................

  if (nonstd) {				# count nonstd symbol entries
     xec[eformat]++;
     if (V) {
        if (N!="")
           print "seq2seq: non-standard nucleotide symbols in \"" name "\"" \
	   	>"/dev/stderr";
        else if (P!="")
           print "seq2seq: non-standard protein symbols in \"" name "\"" \
	   	>"/dev/stderr";
     }
     if (X)
        return 0;
  }

  if (eformat=="how") {			# HOW: seq and ass length check
     if (seqlen!=len) {
        if (V)
           print "seq2seq: HOW sequence length error in \"" name "\"" \
	   	>"/dev/stderr";
        cec[eformat]++;
        return 0;
     }
     else if (asslen!=len) {
        if (V)
           print "seq2seq: HOW assignment length error in \"" name "\"" \
	   	>"/dev/stderr";
        cec[eformat]++;
        return 0;
     }
  }

  # if needed, create empty assignment
  if ((eformat=="fasta") && (OUTPUTFORMAT!="fasta")) {
     for (k=1; k<=len; k++)
	 ass[k]=".";
     asslen = len;    
  }

  # masking: amino acids and/or positions
  if ((M!="") && (asslen==len)) {
     for (k=1; k<=len; k++)
         if ((seq[k] ~ maa) || ((k>=mbottom) && (k<=mtop)))
	    ass[k]="-";
  }

  # truncate
  if (TRUNC>0) {
     if (TRUNC<len)
        len=TRUNC;
  }
  else if (TRUNC<0) {
     if (-1*TRUNC<len) {
        ttn=1;
        for (tt=len+TRUNC+1; tt<=len; tt++) {
	    seq[ttn]=seq[tt];
	    ass[ttn]=ass[tt];
	    ttn++;
	}
        len=-TRUNC;
     }
  }

  # print ...
  oec[OUTPUTFORMAT]++;
  if (OUTPUTFORMAT=="how")
     printhow(OUTDIR,len,name,seq,ass);
  else if (OUTPUTFORMAT=="fasta")
     printfasta(OUTDIR,len,name,seq);
  else if (OUTPUTFORMAT=="tab")
     printtab(OUTDIR,len,name,seq,ass);
  else if (OUTPUTFORMAT=="act")
     printact(OUTDIR,len,name,seq,ass);

  return 0;
}

###############################################################################
function reset() {
###############################################################################

  name		= "";

  len		= 0;
  seqlen	= 0;		# actual sequence length
  asslen	= 0;		# actual assignment length

  corruption	= 0;		# new entry
  nonstd	= 0;		# flag: non-standard symbols detected

  isfastaseq	= 0;		# flag: picking up FASTA sequence
  ishowass	= 0;		# flag: picking up HOW assignment
  ishowseq	= 0;		# flag: picking up HOW sequence

  return 0;

}

###############################################################################
function report() {
###############################################################################

  print "\nseq2seq 1.0 report" >"/dev/stderr";

  # table head
  print "------------------------------------------------------" \
     >"/dev/stderr";
  print "# entries\t\t fasta\t   how\t   tab\t TOTAL" >"/dev/stderr";
  print "------------------------------------------------------" \
     >"/dev/stderr";

  # input
  printf("INPUT\t\t\t%6d\t%6d\t%6d",ec["fasta"],ec["how"],ec["tab"]) \
     >"/dev/stderr";
  printf("\t%6d\n",ec["fasta"]+ec["how"]+ec["tab"]) >"/dev/stderr";

  # corrupted
  printf("Corrupted\t\t%6d\t%6d\t%6d",cec["fasta"],cec["how"],cec["tab"]) \
     >"/dev/stderr";
  printf("\t%6d\n",cec["fasta"]+cec["how"]+cec["tab"]) >"/dev/stderr";

  # zero-length
  printf("Zero-length\t\t%6d\t%6d\t%6d",zec["fasta"],zec["how"],zec["tab"]) \
     >"/dev/stderr";
  printf("\t%6d\n",zec["fasta"]+zec["how"]+zec["tab"]) >"/dev/stderr";

  # too long
  printf("Too long\t\t%6d\t%6d\t%6d",Tec["fasta"],Tec["how"],Tec["tab"]) \
     >"/dev/stderr";
  printf("\t%6d\n",Tec["fasta"]+Tec["how"]+Tec["tab"]) >"/dev/stderr";

  # too short
  printf("Too short\t\t%6d\t%6d\t%6d",Bec["fasta"],Bec["how"],Bec["tab"]) \
     >"/dev/stderr";
  printf("\t%6d\n",Bec["fasta"]+Bec["how"]+Bec["tab"]) >"/dev/stderr";

  # with non-std symbols
  if ((N!="") || (P!="")) {
     printf("With non-std symbols\t%6d\t%6d\t%6d",
            xec["fasta"],xec["how"],xec["tab"]) \
     	>"/dev/stderr";
     printf("\t%6d\n",xec["fasta"]+xec["how"]+xec["tab"]) >"/dev/stderr";
  }

  # with modified name
  printf("With modified name\t%6d\t%6d\t%6d",
         mec["fasta"],mec["how"],mec["tab"]) \
     >"/dev/stderr";
  printf("\t%6d\n",mec["fasta"]+mec["how"]+mec["tab"]) >"/dev/stderr";

  # with truncated name
  if (MAXNAMELEN) {
     printf("With truncated name\t%6d\t%6d\t%6d",
            tec["fasta"],tec["how"],tec["tab"]) \
     	>"/dev/stderr";
     printf("\t%6d\n",tec["fasta"]+tec["how"]+tec["tab"]) >"/dev/stderr";
  }

  # with non-unique name
  printf("With non-unique name\t%6d\t%6d\t%6d",
         dec["fasta"],dec["how"],dec["tab"]) \
     >"/dev/stderr";
  printf("\t%6d\n",dec["fasta"]+dec["how"]+dec["tab"]) >"/dev/stderr";

  # output
  printf("OUTPUT\t\t\t%6d\t%6d\t%6d",oec["fasta"],oec["how"],oec["tab"]) \
     >"/dev/stderr";
  printf("\t%6d\n",oec["fasta"]+oec["how"]+oec["tab"]) >"/dev/stderr";
  print "------------------------------------------------------" \
     >"/dev/stderr";

  return;
}

###############################################################################
function printfasta(OUTDIR,len,name,seq) {
###############################################################################

  # define output file name ...................................................
  if (OUTDIR=="/dev/stdout")
     outfile = OUTDIR;
  else
     outfile = (OUTDIR "/" name);

  # PRINT .....................................................................
  if (PRES)					# header
     print ">" name comment >outfile;
  else
     print ">" name >outfile;

  for (j=1; j<=len; j++) {			# sequence
      printf(seq[j]) >outfile;
      if (!(j%60))
         print "" >outfile;
  }
  if (--j%60)					# '\n' if even multiple of 60
     print "" >outfile;

  # close output file .........................................................
  if (OUTDIR!="/dev/stdout")
     close(outfile);

  return 0;

}

###############################################################################
function printhow(OUTDIR,len,name,seq,ass) {
###############################################################################

  # define output file name ...................................................
  if (OUTDIR=="/dev/stdout")
     outfile = OUTDIR;
  else
     outfile = (OUTDIR "/" name);

  # PRINT .....................................................................
  if (PRES)					# header
     printf("%6d %s %s\n",len,name,comment) >outfile;
  else
     printf("%6d %s\n",len,name) >outfile;

  for (j=1; j<=len; j++) {			# sequence
      printf(seq[j]) >outfile;
      if (!(j%80))
         printf(" %7d\n",j) >outfile;
  }
  if (--j%80)
     print "" >outfile;

  for (j=1; j<=len; j++) {			# assignment
      printf(ass[j]) >outfile;
      if (!(j%80))
         printf(" %7d\n",j) >outfile;
  }
  if (--j%80)
     print "" >outfile;

  # close output file .........................................................
  if (OUTDIR!="/dev/stdout")
     close(outfile);

  return 0;
}

###############################################################################
function printtab(OUTDIR,len,name,seq,ass) {
###############################################################################

  # define output file name ...................................................
  if (OUTDIR=="/dev/stdout")
     outfile = OUTDIR;
  else
     outfile = (OUTDIR "/" name);

  # PRINT .....................................................................
  printf("%s\t",name) >outfile;			# name

  for (j=1; j<=len; j++)			# sequence
      printf(seq[j]) >outfile;

  printf("\t") >outfile;

  for (j=1; j<=len; j++)			# assignment
      printf(ass[j]) >outfile;

  printf("\n") >outfile;

  # close output file .........................................................
  if (OUTDIR!="/dev/stdout")
     close(outfile);

  return 0;
}

###############################################################################
function printact(OUTDIR,len,name,seq,ass) {
###############################################################################

  # define output file name ...................................................
  if (OUTDIR=="/dev/stdout")
     outfile = OUTDIR;
  else
     outfile = (OUTDIR "/" name);

  # PRINT .....................................................................
  for (k=1; k<=len; k++)
      printf("%s %s  %-25s %6d\n",ass[k],seq[k],name,k) >outfile;

  # close output file .........................................................
  if (OUTDIR!="/dev/stdout")
     close(outfile);

  return 0;
}
