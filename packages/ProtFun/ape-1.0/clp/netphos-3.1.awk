#! /usr/freeware/bin/gawk -f

# /usr/cbs/bio/src/ape-1.0/clp/netphos-3.1.awk
# 
# This is NetPhos 3.1 command line parsing script for APE. It is also intended
# as a template script for other servers.
# 
#   #########################################################################
#   ### This script translates the arguments given on the command line to ###
#   ### 'netphos-3.1', a symbolic link to the 'ape' executable,  to a new ###
#   ### argument vector that APE understands.  APE then calls itself with ###
#   ### the new argument vector. 					  ###
#   #########################################################################
#
# Input:	argument vector given to 'netphos-3.1';
# Output:	new argument vector to be given to APE.
# 
# Some options are translated to APE options directly; some are gathered
# in a string to be passed on as "-o str".
# 
# If the words "__VERSION__" or "__HELP__" as given to this script it will
# print the version info or the option summary, respectively.
# 
# VERSION HISTORY
# 	2005 Aug 29	customized for netphos 3.1 by K. Rapacki
# 	2005 Nov 23	'nnhowplayer6' for Linux/i386 replaced by mniel
# 			due to strange bug discovered by Anna L, ucl
# 

BEGIN { # default: APE's own options
	opts = "-m netphos -n netphos-3.1";	# method

	# default: -o option string
	oopts = "-v N=netphos-3.1b";		# name to appear
}

# this will run when -V is given ..............................................
/^__VERSION__/ { print "netphos 3.1b, Nov 2005"; next; }

# this will run when -h is given ..............................................
/^__HELP__/ {

  print "\nnetphos 3.1b, Nov 2005\n";
  
  print "Options:\n";
  print "\t-2\t   run only generic predictions (as in NetPhos 2.0)";
  print "\t-b\t   for each residue display only the best prediction";
  print "\t-c val\t   display only the scores higher than 'val'";
  print "\t-d dir\t   destination directory for graphics files";
  print "\t-g\t   generate GIF graphics";
  print "\t-f gff\t   alternative output format (GFF)";
  print "\t-h\t   print this message";
  print "\t-k\t   run only kinase specific predictions (as in NetPhosK 1.x)";
  print "\t-s,-t,-y   predict only on S, T or Y residues, respectively";
  print "\t-S\t   include the input sequence in the output (GFF only)";
  print "\t-V\t   print version and release date";
  print "\n\t<file>\t   data file(s) to run, FASTA, HOW or TAB\n";

  next;
}

# command line translation ....................................................
{
  for (i=1; i<=NF; i++) {

      if ($i=="-2") {				# -2:	run only NetPhos 2.0
         opts = (opts " -p 'netphos_?'");
	 oopts = (oopts " -v E=1");
      }
      else if ($i=="-b")			# -b:	only best prediction
         oopts = (oopts " -v B=1 -v E=1");
      else if ($i=="-c") {			# -c:	score cut-off
         i++;
         oopts = (oopts " -v C=" $i);
      }
      else if ($i=="-d") {			# -d:	dest. dir for graphics
         i++;
         opts = (opts " -d " $i);
      }
      else if ($i=="-g")			# -g:	generate graphics
         opts = (opts " -g");
      else if ($i=="-f") {			# -f:	output format
         i++;
         oopts = (oopts " -v F=" $i);
      }
      else if ($i=="-h")			# -h:	help!
         opts = (opts " -h");
      else if ($i=="-k")			# -k:	run only NetPhosK 1.x
         opts = (opts " -p " "'netphos_??*'");
      else if ($i=="-s")			# -s:	only serine
         opts = (opts " -i ^S");
      else if ($i=="-t")			# -t:	only threonine
         opts = (opts " -i ^T");
      else if ($i=="-y")			# -y:	only tyrosine
         opts = (opts " -i ^Y");
      else if ($i=="-S")			# -S:	include sequence
         oopts = (oopts " -v S=1");
      else if ($i=="-V")			# -V:	version/date
         opts = (opts " -V");
      else if (index($i,"-")==1)		# -?:	illegal option
         opts = (opts " -h");
      else					# INPUT FILE
         opts = (opts " " $i);
  }

  print "-o \"" oopts "\" " opts;
}

# end of script ===============================================================
