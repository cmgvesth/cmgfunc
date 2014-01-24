#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

#=======================================================================================
#=======================================================================================
#==============================	ERROR HANDLING	========================================
#=======================================================================================
#=======================================================================================

#.............. Print usage ..............
sub usage {	
	print "#==================================================================================\n";
	printf ("%-40s:\t%-50s\n", "# USAGE", "perl CMGfunc_createBLAST2GOannot.pl -res <CMGfunc res file> -cut <comparison score>");	
	printf ("%-40s:\t%-50s\n", "# USAGE", "-cut is a number between 0 and 1, where 1 would be a perfect match to the fucntional model, default is 0.8");	
	printf ("%-40s:\t%-50s\n", "# EXAMPLE", "perl CMGfunc_createBLAST2GOannot.pl -res CAFA.proteins.fsa.tab.vector.res");
	print "#==================================================================================\n";
	exit( 1 );	}

#.............. Inputs ..............
my ($res, $h, $scorecut) = (undef,undef, 0);
my $INFO_path		= "/home/cmgfunc/CMGfunc/INFO";

#.............. Get options ..............
&GetOptions ("res:s" =>  \$res, "cut=s" =>  \$scorecut, "h|?" => \$h);
if (defined $h) { &usage }

unless (defined $res)	{	
	print "#================================================================================================================================\n";	
	printf ("%-40s:\t%-50s\n", "# ERROR", "CMGfunc RES file not defined, example: perl CMGfunc_createBLAST2GOannot.pl -res file.fsa.tab.vector.res");	
	&usage;	exit; }

unless ($scorecut >=0 and $scorecut <= 1) {
	print "#================================================================================================================================\n";	
	printf ("%-40s:\t%-50s\n", "# ERROR", "Cutoff score is not a number between 0 and 1: $scorecut");	
	&usage;	exit; }

unless (&checkFileTypeRESULTS($res) eq "TRUE") {	
	print "#================================================================================================================================\n";	
	printf ("%-40s\t==>\t%-50s\n", "# ERROR", "File, $res, is not CMGfunc RES format");	
	&usage;	exit }

#=======================================================================================
#=======================================================================================
#==================================	INFORMATION	========================================
#=======================================================================================
#=======================================================================================

#.............. Print information to user ..............
my $res_seq_count	= `awk 'BEGIN{FS="\t"}{print \$2}' $res | uniq | wc | awk '{print \$1}'`; chomp $res_seq_count;
print "#================================================================================================================================\n";
printf ("%-40s:\t%-80s|\n", '# CMGfunc RES file', 	$res);
printf ("%-40s:\t%-80s|\n", '# Seqs. in RES file',	$res_seq_count); 
print "#================================================================================================================================\n";

#=======================================================================================
#=======================================================================================
#=======================================	MAIN	================================
#=======================================================================================
#=======================================================================================

my ($net, $seq);
my (%GOHash, %freqHash);

printf ("%-40s\t==>\t%-50s\n", "# RUNNING createBLAST2GOannot", "will print to screen");

open(FHres,$res);
while( my $line = <FHres> )  {
    chomp $line;
	my (%net_go);
	 
	next if $line =~ m/^#/ or $line =~ m/group_concat/;
	my @elements = split("\t", $line);
	my ($seq, $net, $score) = ($elements[1], $elements[3], $elements[5]);
	next if $score < $scorecut;

	if ($net =~ m/\_clan\_/) {
		my ($net_nr, $go_type, $type_go, $go_terms) = ('','NA','NA');
		my (@pfams, @go_tmp);

		$net_nr = (split("\\.", $net))[0];
		$net_nr =~ s/proteindata\_clan\_//g;
		$go_terms = `cut -f 3,6 /home/cmgfunc/CMGfunc/INFO/info_pfam_clannr_clandesc_gos_godescs.txt | grep -P "^$net_nr\t" | cut -f 2 | perl -pe 's/::/\n/g' | sort -u | perl -pe 's/\n/::/g'`;

		@go_tmp = split("::", $go_terms);
		#print "$net_nr :: ",scalar(@go_tmp),"\n";
		
		if ($go_terms ne "NA::") {
			for my $go (@go_tmp) {
				next if ($go eq "NA");
				#if ($go eq "NA") { print "$seq\t$net_nr\tNA\t$score\n"; next }
				#$go_type = `grep -P $go /home/cmgfunc/CMGfunc/INFO/info_goid_desc_type_go.txt | cut -f 3`; chomp $go_type;
				#if		($go_type =~ m/molecular/)	{	$type_go = $go."_MF" } 
				#elsif	($go_type =~ m/cellular/)	{	$type_go = $go."_CC" } 
				#elsif	($go_type =~ m/biological/) {	$type_go = $go."_BP" }
				$type_go = "GO:".$go;
				print "$seq\t$type_go\n";
			}	
		}	

	}
	elsif ($net =~ m/\_arch\_/) {
		my ($net_nr, $pfam, $go_type, $type_go, $go_terms) = ('','','NA','NA');
		my (@pfams, @go_tmp);

		$net_nr = (split("\\.", $net))[0];
		$net_nr =~ s/proteindata\_arch\_//g;
		@pfams = split("\_", $net_nr);

		for $pfam (@pfams) {
			$go_terms = `grep -P -m 1 "^$pfam\t" /home/cmgfunc/CMGfunc/INFO/info_pfam_clannr_clandesc_gos_godescs.txt | cut -f 6`; chomp $go_terms;
			@go_tmp = split("::", $go_terms);
			if ($go_terms ne "NA") {
				for my $go (@go_tmp) {
					#if ($go eq "NA") { print "$seq\t$net_nr\tNA\t$score\n"; next }
					#$go_type = `grep -P $go /home/cmgfunc/CMGfunc/INFO/info_goid_desc_type_go.txt | cut -f 3`; chomp $go_type;
					#if		($go_type =~ m/molecular/)	{	$type_go = $go."_MF" } 
					#elsif	($go_type =~ m/cellular/)	{	$type_go = $go."_CC" } 
					#elsif	($go_type =~ m/biological/) {	$type_go = $go."_BP" }
					$type_go = "GO:".$go;
					print "$seq\t$type_go\n";
				}	
			}	
		}
	}
}


#=======================================================================================
#=======================================================================================
#=======================================	SUBS	================================
#=======================================================================================
#=======================================================================================

#################################################
sub checkFileTypeRESULTS {
#################################################
	my $status			= "TRUE";
	my $res_filename	= $_[0];
	my $i				= 0;

	# Open the TAB file and test each vector for problems
	open(FP0,$res_filename);
	while( my $line = <FP0> )  {
		next if $line =~ /^#/; 
		my @elements = split("\t", $line);

		unless (scalar(@elements) == 6)	{ 
			printf ("%-40s\t==>\t%-50s%-10s\n", "# ERROR checkFileTypeRESULTS", "RES vectore, $res_filename, line $i does not have the right number of elements");
			#print scalar(@elements), "\n";
			$status = "FALSE"; return $status; }

		foreach my $item (@elements)		{ 
			if ($item eq '') { printf ("%-40s\t==>\t%-50s\n", "# ERROR checkFileTypeRESULTS", "RES vectore, $res_filename, in line $i has empty elements");		
			$status = "FALSE"; return $status; }}

		$i++;
	}
	close(FP0);
	return $status;
}
