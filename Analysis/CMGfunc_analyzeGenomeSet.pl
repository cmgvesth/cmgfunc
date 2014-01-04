#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Data::Dumper;
use Bio::SeqIO;
=put
#############################################################################
# Description
#############################################################################

#--------------------------------------------------------------
# Input source
#--------------------------------------------------------------
CMGfunc_testNetworks.py -i file.fsa.tab.vec -n /home/cmgfunc/CMGfunc/NETS/ -m 0.8

OR

perl CMGfunc.pl -fasta file.fsa -net /home/cmgfunc/CMGfunc/NETS -mse 0.8

#--------------------------------------------------------------
# Input
#--------------------------------------------------------------

head data_seq_trial_clan.fsa.tab.vec.res
Testseq:	10002486_PF02686_NULL 	Network:	proteindata_clan_157.tab.A.net 	Score:	0.87336622785
Testseq:	10002486_PF02686_NULL 	Network:	proteindata_clan_190.tab.A.net 	Score:	1.05753914864
Testseq:	10002486_PF02686_NULL 	Network:	proteindata_clan_214.tab.A.net 	Score:	0.85312832223
Testseq:	10002486_PF02686_NULL 	Network:	proteindata_clan_255.tab.A.net 	Score:	0.858665730197
Testseq:	10002486_PF02686_NULL 	Network:	proteindata_clan_509.tab.A.net 	Score:	0.983983750904

#--------------------------------------------------------------
# Output
#--------------------------------------------------------------
Freq. table and plots

#--------------------------------------------------------------
# USAGE
#--------------------------------------------------------------

perl CMGfunc_analyzeGenomeSet.pl -res file.fsa.tab.vec.res
=cut

#=======================================================================================
#=======================================================================================
#==============================	ERROR HANDLING	========================================
#=======================================================================================
#=======================================================================================

#.............. Inputs ..............
my ($res, $cut, $h);
my $CMGfunc_path = `locate CMGfunc | grep -P "CMGfunc/indirect" | grep -P -v "CMGfunc/indirect." | grep -v MySQL`; chomp $CMGfunc_path;

unless (-d $CMGfunc_path)	{	#.............. CMGfunc path ..............
	printf ("%-40s:\t%-50s\n", "# ERROR checkPATHS", "nothing called $CMGfunc_path or is not a directory"); 
	printf ("%-40s:\t%-50s\n", "# ERROR checkPATHS", "CMGfunc_directory : $CMGfunc_path");		
	exit }
#.............. Get options ..............
&GetOptions ("res:s" =>  \$res, "cut:f" => \$cut, "h|?" => \$h);

if (defined $h) { &usage }

#.............. Input files and paths exists ..............
unless (defined $res)	{	
	printf ("%-40s:\t%-50s\n", "# ERROR", "CMGfunc RES file not defined, example: perl CMGfunc_analyzeGenome.pl -res genomeA.proteins.fsa.tab.vector.res");	&usage;	exit; }

unless (defined $cut)	{
	printf ("%-40s:\t%-50s\n", "# WARNING", "Frequency cutoff not defined (-cut), default used 10");	
	$cut	= 10;	}


#.............. Print usage ..............
sub usage {	
	#print "#==================================================================================\n";
	printf ("%-40s:\t%-50s\n", "# USAGE", "perl CMGfunc_analyzeGenome.pl -res <CMGfunc res file> -mse <MSE value>");	
	printf ("%-40s:\t%-50s\n", "# EXAMPLE", "perl CMGfunc_analyzeGenome.pl -res genomeA.proteins.fsa.tab.vector.res -cut 20");
	printf ("%-40s:\t%-50s\n", "# INFO", "Frequency cutoff is optional, default is 10, cutoff is used to filter network comparison results");
	#print "#==================================================================================\n";
	exit( 1 );	}

#.............. Test that paths exists or exit ..............
&checkPATHS();

#.............. Test that file is CMGfunc res file or exit ..............
unless (&checkFileTypeRESULTS($res) eq "TRUE")		{	printf ("%-40s\t==>\t%-50s\n", "# ERROR", "File, $res, is not CMGfunc RES format");	exit }

#=======================================================================================
#=======================================================================================
#==================================	INFORMATION	========================================
#=======================================================================================
#=======================================================================================
my $res_seq_count	= `awk '{print \$1}' $res | uniq | wc | awk '{print \$1}'`; chomp $res_seq_count;

#.............. Print information to user ..............
print "#================================================================================================================================\n";
printf ("%-40s:\t%-80s|\n", '# CMGfunc RES file', 			$res);
printf ("%-40s:\t%-80s|\n", '# Seqs. in RES file', 			$res_seq_count); 
printf ("%-40s:\t%-80s|\n", '# Frequency cutoff',			$cut);
print "#================================================================================================================================\n";

#=======================================================================================
#=======================================================================================
#=======================================	MAIN	================================
#=======================================================================================
#=======================================================================================

#if ($clean == 0) {
#	printf ("%-40s\t==>\t%-50s\n", "# OPTION -clean defined", "REMOVING files $res.tab.vector.res.func");
#	`rm -f $res.tab.vector.res.func`;
#}

&analyzeResults($res);	#.............. Analyze results ..............
unless (&checkFileTypeFUNCTIONS($res.".func") eq "TRUE")		{ exit }

&plotResults($res.".func");

#=======================================================================================
#=======================================================================================
#==========================	CALCULATION SUBRUTINES	====================================
#=======================================================================================
#=======================================================================================
#.............. analyzeResults ..............
sub analyzeResults {
	my $res_filename = $_[0];	
	my (%Hres, $net, $seq, $score);

	printf ("%-40s\t==>\t%-50s\n", "# RUNNING analyzeResults", "will create file: $res_filename.func");

	open(FHres,$res_filename);
	while( my $line = <FHres> )  {
		my @elements = split("\t", $line);
		$seq = $elements[1];
		$net = $elements[3];
		$score = $elements[5];
		$Hres{$net} ++;
	}
	
	while (my ($n, $s) = each(%Hres)) {
		print "$n :: $s\n" if $s >= $cut;	
	}

	unless (-e $cd."/".$res.".func")	{ 
		printf ("%-40s\t==>\t%-50s\n", "# ERROR analyzeResults", "$CMGfunc_path/CMGfunc_analyzeResults.pl did not return file $res.func");	exit }
}

#.............. analyzeResults ..............
sub plotResults {
	my $res_filename = $_[0];	
	printf ("%-40s\t==>\t%-50s\n", "# RUNNING plotResults", "will create file: $res_filename.func.pdf");


	unless (-e $res_filename.".func.pdf")	{ 
		printf ("%-40s\t==>\t%-50s\n", "# ERROR plotResults", "$CMGfunc_path/CMGfunc_analyzeResults.pl did not return file $res_filename.func.pdf");	exit }
}



#=======================================================================================
#=======================================================================================
#===========================	FILETYPE SUBRUTINES	====================================
#=======================================================================================
#=======================================================================================

#.............. checkFileTypeFUNCTIONS ..............
sub checkFileTypeFUNCTIONS {
	my $func_filename	= $_[0];
	my $status = "TRUE";
	return $status;
}

#.............. checkFileTypeRESULTS ..............
sub checkFileTypeRESULTS {
	my $status			= "TRUE";
	my $res_filename	= $_[1];
	my $fsa_filename	= $_[0];
	my $net_library		= $_[2];
	my $i				= 0;
	#my $fsa_seq_count	= `grep -c '>' $fsa_filename`;
	#my $res_line_count	= `wc $res_filename | awk '{print \$1}'`;
	#my $net_files		= `ls -1 $NET_path | wc | awk '{print \$1}'`;

	# Open the TAB file and test each vector for problems
	open(FP0,$res_filename);
	while( my $line = <FP0> )  {
		my @elements = split("\t", $line);

		unless (scalar(@elements) == 6)	{ 
			printf ("%-40s\t==>\t%-50s\n", "# ERROR checkFileTypeRESULTS", "RES vectore, $res_filename, in line $i does not have the right number of elements");
			$status = "FALSE"; return $status; }

		foreach my $item (@elements)		{ 
			if ($item eq '') { printf ("%-40s\t==>\t%-50s\n", "# ERROR checkFileTypeRESULTS", "RES vectore, $res_filename, in line $i has empty elements");		
			$status = "FALSE"; return $status; }}

		#unless ($elements[1] =~ /[a-zA-Z]/)	{ 
		#	printf ("%-40s\t==>\t%-50s\n", "# ERROR checkFileTypeRESULTS", "RES vectore, $res_filename, first element is a number not a seq. name");			
		#	$status = "FALSE"; return $status; }
		$i++;
	}
	close(FP0);
	return $status;
}
#=======================================================================================
#=======================================================================================
#===========================	PATH AND FILE CHECK SUBS	============================
#=======================================================================================
#=======================================================================================

#.............. checkPATHS ..............
sub checkPATHS {
	unless (-d $CMGfunc_path)	{	#.............. CMGfunc path ..............
		printf ("%-40s:\t%-50s\n", "# ERROR checkPATHS", "nothing called $CMGfunc_path or is not a directory"); 
		printf ("%-40s:\t%-50s\n", "# ERROR checkPATHS", "CMGfunc_directory : $CMGfunc_path");		
		exit }
}
