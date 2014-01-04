#! /usr/bin/nawk -f

# 
# loadenv.awk
#
# This script prints environment setting commands based on a configuration file
# given as argument. The file has to have the following syntax:
# 
# 	<varname>	<value>
# 
# The values may contain whitespace.
# 
# SYNOPSIS
# 	eval `<full_path>/loadenv cffile`
# 
# VERSIONS		2004 Sep 14	launch, K. Rapacki
# 


{ gsub("#.*$","",$0); }				# strip comments

! NF { next; }					# skip empty lines

{
  varname = toupper($1);			# get variable name
  gsub("[^A-Z]","_",varname);			# make sure it is healthy ...
  sub("[^ \t]*[ \t]*","",$0);			# get value
  sub("[ \t]*$","",$0);				# strip trailing whitespace

  print "setenv " varname " \"" $0 "\";";	# print setenv command
}

# end of script ===============================================================
