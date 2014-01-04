#!/usr/bin/perl

use strict;
use warnings;
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
# Comparison started for file:  proteindata_clan_101.tab.A.test.tab.vec
Testseq:	name1 _nr_ 1 	Network:	proteindata_clan_65.tab.C.net 	Score:	0.833411100928
Testseq:	name1 _nr_ 1 	Network:	proteindata_clan_122.tab.A.mse_0.0099.net 	Score:	0.907742353911
Testseq:	name1 _nr_ 1 	Network:	proteindata_clan_101.tab.A.mse_0.0121.net 	Score:	0.995877724143
Testseq:	name1 _nr_ 1 	Network:	proteindata_clan_101.tab.C.net 	Score:	1.06325840171
Testseq:	name1 _nr_ 3 	Network:	proteindata_clan_65.tab.C.net 	Score:	0.806739932733
Testseq:	name1 _nr_ 3 	Network:	proteindata_clan_101.tab.A.mse_0.0121.net 	Score:	1.00185420067

#--------------------------------------------------------------
# Output
#--------------------------------------------------------------
Freq. table and plots

#--------------------------------------------------------------
# USAGE
#--------------------------------------------------------------
perl CMGfunc_analyzeGenome.pl -res file.fsa.tab.vec.res

cmgfunc@cmgfunc-VirtualBox:~/CMGfunc$ for x in RESTrials/*res
do
perl CMGfunc_analyzeGenome.pl -res $x 
echo $x
done

=cut

#=======================================================================================
#=======================================================================================
#==============================	ERROR HANDLING	========================================
#=======================================================================================
#=======================================================================================

#.............. Inputs ..............
my ($res, $cut, $h);
my $CMGfunc_path = `locate CMGfunc | grep -P "CMGfunc/indirect" | grep -P -v "CMGfunc/indirect." | grep -v MySQL`; chomp $CMGfunc_path;
my $INFO_path = `locate CMGfunc | grep -P "CMGfunc/indirect" | grep -P -v "CMGfunc/indirect." | grep -v MySQL`; chomp $CMGfunc_path;

unless (-d $CMGfunc_path)	{	#.............. CMGfunc path ..............
	printf ("%-40s:\t%-50s\n", "# ERROR checkPATHS", "nothing called $CMGfunc_path or is not a directory"); 
	printf ("%-40s:\t%-50s\n", "# ERROR checkPATHS", "CMGfunc_directory : $CMGfunc_path");		
	exit }

#.............. Get options ..............
&GetOptions ("res:s" =>  \$res, "cut:f" => \$cut, "h|?" => \$h);
if (defined $h) { &usage }

#.............. Input files and paths exists ..............
unless (defined $res)	{	
	printf ("%-40s:\t%-50s\n", "# ERROR", "CMGfunc RES file not defined, example: perl CMGfunc_analyzeGenome.pl -res file.fsa.tab.vector.res");	&usage;	exit; }

unless (defined $cut)	{
	printf ("%-40s:\t%-50s\n", "# WARNING", "Frequency cutoff not defined (-cut), default used 10");	
	$cut	= 10;	}


#.............. Print usage ..............
sub usage {	
	#print "#==================================================================================\n";
	printf ("%-40s:\t%-50s\n", "# USAGE", "perl CMGfunc_analyzeGenome.pl -res <CMGfunc res file> -cut <MSE value>");	
	printf ("%-40s:\t%-50s\n", "# EXAMPLE", "perl CMGfunc_analyzeGenome.pl -res file.proteins.fsa.tab.vector.res -cut 20");
	printf ("%-40s:\t%-50s\n", "# INFO", "Frequency cutoff is optional, default is 10, cutoff is used to desplay function distribution results");
	#print "#==================================================================================\n";
	exit( 1 );	}

#.............. Test that file is CMGfunc res file or exit ..............
unless (&checkFileTypeRESULTS($res) eq "TRUE")		{	printf ("%-40s\t==>\t%-50s\n", "# ERROR", "File, $res, is not CMGfunc RES format");	exit }

#=======================================================================================
#=======================================================================================
#==================================	INFORMATION	========================================
#=======================================================================================
#=======================================================================================

#.............. Print information to user ..............
my $res_seq_count	= `awk 'BEGIN{FS="\t"}{print \$2}' $res | uniq | wc | awk '{print \$1}'`; chomp $res_seq_count;print "#================================================================================================================================\n";
printf ("%-40s:\t%-80s|\n", '# CMGfunc RES file', 	$res);
printf ("%-40s:\t%-80s|\n", '# Seqs. in RES file',	$res_seq_count); 
printf ("%-40s:\t%-80s|\n", '# Frequency cutoff',	$cut);
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
#unless (&checkFileTypeFUNCTIONS($res.".func") eq "TRUE")		{ exit }

#&plotResults($res.".func");

#=======================================================================================
#=======================================================================================
#=======================================	SUBRUTINES	================================
#=======================================================================================
#=======================================================================================

#################################################
sub analyzeResults {
#################################################
	my $res_filename = $_[0];	
	my ($net, $seq);
	my (%GOHash, %lowGOHash) = ((),());
	my (%ClanHash, %lowClanHash, %ArchClanHash, %lowArchClanHash, %ArchNoClanHash, %lowArchNoClanHash);

	printf ("%-40s\t==>\t%-50s\n", "# RUNNING analyzeResults", "will create file: $res_filename.func");

	open(FHres,$res_filename);
	while( my $line = <FHres> )  {
        chomp $line;
		next if $line =~ m/^#/ or $line =~ m/group_concat/;
		my @elements = split("\t", $line);

		my ($seq, $net, $score) = ($elements[1], $elements[3], $elements[5]);
		
		if ($net =~ m/\_clan\_/) { 
			#print "PRINT-- ",Dumper(%ClanHash), " --\n";
			#print "$net :: $score :: ", %GOHash ," :: ",%lowGOHash, "\n";
			&clanNetworks($net, $score, \%GOHash, \%lowGOHash, \%ClanHash, \%lowClanHash);
		}
		elsif ($net =~ m/\_arch\_/) {
			&archClanNetworks($net, $score, \%GOHash, \%lowGOHash, \%ArchClanHash, \%lowArchClanHash);
		}
		elsif ($net =~ m/\_archNOclan\_/) {
			&archNoClanNetworks($net, $score, \%GOHash, \%lowGOHash, \%ArchNoClanHash, \%lowArchNoClanHash);
		}

	} # WHILE RESLINE END  ---------------------------------------------
	#print "PRINT-- ",Dumper(%ArchNoClanHash), " --\n";

	#----------------------------------------------------------------
	# PRINTS	-----------------------------------------------------
	#----------------------------------------------------------------

	open(FHtable, ">", $res_filename.".table");
	open(FHplottable, ">", $res_filename.".plottable");
	
	&print_clanNetworks(\%ClanHash, "high");
	&print_archClanNetworks(\%ArchClanHash, "high");
	&print_archNoClanNetworks(\%ArchNoClanHash, "high");
	&print_GOs(\%GOHash, "high");

	&print_clanNetworks(\%lowClanHash, "low");
	&print_archClanNetworks(\%lowArchClanHash, "low");
	&print_archNoClanNetworks(\%lowArchNoClanHash, "low");
	&print_GOs(\%lowGOHash, "low");

	close FHres;
	close FHtable;
	close FHplottable;

}	# SUB RUTINE END  ---------------------------------------------

#=======================================================================================
#=======================================================================================
#==========================	CALCULATION SUBRUTINES	====================================
#=======================================================================================
#=======================================================================================

#################################################
sub clanNetworks {
#################################################
	my ($net, $score, $GOHash, $lowGOHash, $ClanHash, $lowClanHash) = @_;

	# HASH of clan networks	---------------------------------------------
	my $netnr = (split("_",(split("\\.", $net))[0]))[-1];	
	my %go_per_arch;

	# GO freq. counts ---------------------------------------------
	my $go_terms = `grep -P '\t$netnr\t' /home/cmgfunc/CMGfunc/INFO/info_gos_clan_pfam.txt | cut -f 1`; chomp $go_terms;
	my @gos = split(",", $go_terms);
	for my $g (0 .. $#gos) {
		my $go_strip = $gos[$g]; $go_strip =~ s/GO://g;	
		$$GOHash{$go_strip}++ if ($score >= 0.8);
		$$lowGOHash{$go_strip}++ if ($score < 0.8);
		$go_per_arch{$go_strip}++ 
	}
	# Seperate high and low scoring results --------------------------------
	if		(scalar(keys %go_per_arch) > 0 and $score >= 0.8)	{ $$ClanHash{$netnr}{scalar(keys %go_per_arch)} ++ }
	elsif	(scalar(keys %go_per_arch) > 0 and $score < 0.8)	{ $$lowClanHash{$netnr}{scalar(keys %go_per_arch)} ++ }
	elsif	(scalar(keys %go_per_arch) <= 0 and $score < 0.8)	{ $$lowClanHash{$netnr}{"0"} ++ }
	else	{ $$ClanHash{$netnr}{"0"} ++ }

	#return ($GOHash, $lowGOHash, $ClanHash, $lowClanHash);
}
#################################################
#################################################

#################################################
sub archClanNetworks {
#################################################
	my ($net, $score, $GOHash, $lowGOHash,$ArchClanHash, $lowArchClanHash) = @_;

	# HASH of Archclan networks	-------------------------------------
	my $netarch = (split("\\.", $net))[0]; $netarch =~ s/proteindata_arch_//g;
	my @pfams = (split("_", $netarch));
	my %go_per_arch;

	# GO freq. counts ---------------------------------------------
	for my $p (@pfams) {		
		my @go_terms = `grep -P '^$p\t' /home/cmgfunc/CMGfunc/INFO/info_go_pfam.txt | cut -f 2`; chomp @go_terms;
		for my $g (0 .. $#go_terms) {	
			my $go_strip = $go_terms[$g]; $go_strip =~ s/GO://g;	
			$$GOHash{$go_strip}++ if ($score >= 0.8);
			$$lowGOHash{$go_strip}++ if ($score < 0.8);
			$go_per_arch{$go_strip}++ 
		}
	}
	# Seperate high and low scoring results --------------------------------
	if		(scalar(keys %go_per_arch) > 0 and $score >= 0.8)	{ $$ArchClanHash{$netarch}{scalar(keys %go_per_arch)} ++ }
	elsif	(scalar(keys %go_per_arch) > 0 and $score < 0.8)	{ $$lowArchClanHash{$netarch}{scalar(keys %go_per_arch)} ++ }
	elsif	(scalar(keys %go_per_arch) <= 0 and $score < 0.8)	{ $$lowArchClanHash{$netarch}{"0"} ++ }
	else	{ $$ArchClanHash{$netarch}{"0"} ++ }

}
#################################################
#################################################

#################################################
sub archNoClanNetworks {
#################################################
	my ($net, $score, $GOHash, $lowGOHash, $ArchNoClanHash, $lowArchNoClanHash) = @_;

	# HASH of ArchNoclan networks	---------------------------------
	my $netarch = (split("\\.", $net))[0]; $netarch =~ s/proteindata_archNOclan_//g;
	my @pfams = (split("_", $netarch));
	my %go_per_arch;

	# GO freq. counts ---------------------------------------------
	for my $p (@pfams) {		
		my @go_terms = `grep -P '^$p\t' /home/cmgfunc/CMGfunc/INFO/info_go_pfam.txt | cut -f 2`; chomp @go_terms;
		for my $g (0 .. $#go_terms) {	
			my $go_strip = $go_terms[$g]; $go_strip =~ s/GO://g;	
			$$GOHash{$go_strip}++ if ($score >= 0.8); 
			$$lowGOHash{$go_strip}++ if ($score < 0.8); 
			$go_per_arch{$go_strip}++ 
		}
	}
	if		(scalar(keys %go_per_arch) > 0 and $score >= 0.8)	{ $$ArchNoClanHash{$netarch}{scalar(keys %go_per_arch)} ++ }
	elsif	(scalar(keys %go_per_arch) > 0 and $score < 0.8)	{ $$lowArchNoClanHash{$netarch}{scalar(keys %go_per_arch)} ++ }
	elsif	(scalar(keys %go_per_arch) <= 0 and $score < 0.8)	{ $$lowArchNoClanHash{$netarch}{"0"} ++ }
	else	{ $$ArchNoClanHash{$netarch}{"0"} ++ }
}
#################################################
#################################################

#=======================================================================================
#=======================================================================================
#==========================	PRINTING SUBRUTINES	====================================
#=======================================================================================
#=======================================================================================

#################################################
sub print_clanNetworks {
#################################################
	my ($ClanHash, $type) = @_;

	print  				"#TYPE\t$type:CLAN\n";
	print FHplottable	"#TYPE\t$type:CLAN\tCLAN name\tFreq\tNr. GO\n";
	print FHtable		"#TYPE\t$type:CLAN\tCLAN name\tFreq\tNr. GO\tClanDesc\n";

	foreach my $n (keys %$ClanHash)	{	foreach my $d (keys %{$$ClanHash{$n}})	{	
		my $clan_name = `grep -P '^$n\t' /home/cmgfunc/CMGfunc/INFO/info_clan_desc_pfam_com.txt | cut -f 2`; chomp $clan_name;
		my $clan_desc = `grep -P '^$n\t' /home/cmgfunc/CMGfunc/INFO/info_clan_desc_pfam_com.txt | cut -f 5`; chomp $clan_desc;

		print FHtable		"$type:CLAN $n\t$clan_name\t$$ClanHash{$n}{$d}\t$clan_desc\t$d\n";
		print FHplottable	"$type:CLAN\t$n\t$clan_name\t$$ClanHash{$n}{$d}\t$d\n";
	}} 
}
#################################################
#################################################

#################################################
sub print_archClanNetworks {
#################################################
	my ($ArchClanHash, $type) = @_;

	print  				"#TYPE\t$type:ARCH\n";
	print FHtable		"#TYPE\t$type:ARCH\tCLAN name\tFreq\tNr. GO\tClanDesc\tPfamDescs\tGODescs\tGOs\n";
	print FHplottable	"#TYPE\t$type:ARCH\tCLAN name\tFreq\tNr. GO\n";

	foreach my $n ( keys %$ArchClanHash )	{	foreach my $d ( keys %{$$ArchClanHash{$n}} )	{	
		my (@pfam_go_types, @pfam_go_descs, @pfam_descs, @pfam_gos) = (),(),(),();
	    my $clannr 		= `cut -f 1,3 /home/cmgfunc/CMGfunc/INFO/info_clan_pfams_uarch_nrpfam_nrarchs.txt | grep -P ',$n,' |  cut -f 1`; chomp $clannr;
		my $clan_name	= `grep -P '^$clannr\t' /home/cmgfunc/CMGfunc/INFO/info_clan_desc_pfam_com.txt | cut -f 2`; chomp $clan_name;
		my $clan_desc	= `grep -P '^$clannr\t' /home/cmgfunc/CMGfunc/INFO/info_clan_desc_pfam_com.txt | cut -f 5`; chomp $clan_desc;
		my @pfams		= (split("_", $n)); 

		for my $p (@pfams) {		
			my $pfam_desc = `grep -P '^$p\t' /home/cmgfunc/CMGfunc/INFO/info_pfamdomain_desc.txt | cut -f 2`; chomp $pfam_desc;
			push(@pfam_descs, $pfam_desc);
			@pfam_gos = `grep -P '^$p\t' /home/cmgfunc/CMGfunc/INFO/info_go_pfam.txt | cut -f 2`; chomp @pfam_gos;
			
			if (scalar(@pfam_gos) < 1) { @pfam_go_descs = ("none") ; @pfam_gos = ("none"); @pfam_go_types = ("none"); }
	
			else { for my $g (@pfam_gos) {
				my $pfam_go_desc	= `grep -P '$g' /home/cmgfunc/CMGfunc/INFO/info_goid_desc_type_go.txt | cut -f 2`;	chomp $pfam_go_desc; $pfam_go_desc =~ s/^ //g;
				my $pfam_go_type	= `grep -P '$g' /home/cmgfunc/CMGfunc/INFO/info_goid_desc_type_go.txt | cut -f 3`;	chomp $pfam_go_type; $pfam_go_type =~ s/^ //g;
				push(@pfam_go_descs, $pfam_go_desc) unless $pfam_go_desc eq '';				
				push(@pfam_go_types, $pfam_go_type) unless $pfam_go_type eq '';				
			}}

			print FHtable		"$type:ARCHCLAN $n\t$clan_name\t$$ArchClanHash{$n}{$d}\t$d\t$clan_desc\t", join(";",@pfam_descs),"\t", join(";",@pfam_go_descs),"\t", join(";",@pfam_gos),"\n";
			print FHplottable	"$type:ARCHCLAN $n\t$clan_name\t$$ArchClanHash{$n}{$d}\t$d\n";				
		}  
	}} 
}
#################################################
#################################################

#################################################
sub print_archNoClanNetworks {
#################################################
	my ($ArchNoClanHash, $type) = @_;

	print  				"#TYPE\t$type:ARCHnoClan\n";
	print FHtable		"#TYPE\t$type:ARCH\tFreq\tNr. GO\tPfamDescs\tGODescs\tGOs\n";
	print FHplottable	"#TYPE\t$type:ARCH\tFreq\tNr. GO\n";

	foreach my $n ( keys %$ArchNoClanHash )	{	foreach my $d ( keys %{$$ArchNoClanHash{$n}} )	{	

		my @pfams = (split("_", $n)); 
		my (@pfam_go_types, @pfam_go_descs, @pfam_descs, @pfam_gos);

		for my $p (@pfams) {		
			my $pfam_desc = `grep -P '^$p\t' /home/cmgfunc/CMGfunc/INFO/info_pfamdomain_desc.txt | cut -f 2`; chomp $pfam_desc;
			push(@pfam_descs, $pfam_desc);
			@pfam_gos = `grep -P '^$p\t' /home/cmgfunc/CMGfunc/INFO/info_go_pfam.txt | cut -f 2`; chomp @pfam_gos;
			
			if (scalar(@pfam_gos) < 1) { @pfam_go_descs = ("none") ; @pfam_gos = ("none"); @pfam_go_types = ("none"); }
	
			else { for my $g (@pfam_gos) {
				my $pfam_go_desc	= `grep -P '$g' /home/cmgfunc/CMGfunc/INFO/info_goid_desc_type_go.txt | cut -f 2`;	chomp $pfam_go_desc; $pfam_go_desc =~ s/^ //g;
				my $pfam_go_type	= `grep -P '$g' /home/cmgfunc/CMGfunc/INFO/info_goid_desc_type_go.txt | cut -f 3`;	chomp $pfam_go_type; $pfam_go_type =~ s/^ //g;
				push(@pfam_go_descs, $pfam_go_desc) unless $pfam_go_desc eq '';				
				push(@pfam_go_types, $pfam_go_type) unless $pfam_go_type eq '';				
			}}
			print FHtable		"$type:ARCHnoCLAN $n\t$$ArchNoClanHash{$n}{$d}\t$d\t", join(";",@pfam_descs),"\t", join(";",@pfam_go_descs),"\t", join(";",@pfam_gos),"\n";
			print FHplottable	"$type:ARCHnoCLAN $n\t$$ArchNoClanHash{$n}{$d}\t$d\n";
		}
	}} 
}
#################################################
#################################################

#################################################
sub print_GOs {
#################################################
	my ($GOHash, $type) = @_;

	print  "#TYPE\t$type:GO\n";
	foreach my $g ( keys %$GOHash )	{
		my $go_desc			= `grep -P '$g' /home/cmgfunc/CMGfunc/INFO/info_goid_desc_type_go.txt | cut -f 2`;	chomp $go_desc; $go_desc =~ s/^ //g;
		my $go_id 			= `grep -P '$g' /home/cmgfunc/CMGfunc/INFO/info_goid_desc_type_go.txt | cut -f 1`;	chomp $go_id;
		my $go_long_desc	= `grep -P '^$go_id\t' /home/cmgfunc/CMGfunc/INFO/info_clanid_desc.txt | cut -f 2`; chomp $go_long_desc;
		my $go_type 		= `grep -P '$g' /home/cmgfunc/CMGfunc/INFO/info_goid_desc_type_go.txt | cut -f 3`;	chomp $go_type; $go_type =~ s/^ //g;

		print FHtable		"$type:GO\t$g\t$go_desc\t$go_type\t$$GOHash{$g}\t$go_long_desc\n";
		print FHplottable	"$type:GO\t$g\t$go_desc\t$go_type\t$$GOHash{$g}\n";
	}
}
#################################################
#################################################


#################################################
sub plotResults {
#################################################
	my $res_filename = $_[0];	
	printf ("%-40s\t==>\t%-50s\n", "# RUNNING plotResults", "will create file: $res_filename.func.pdf");
	unless (-e $res_filename.".func.pdf")	{ 
		printf ("%-40s\t==>\t%-50s\n", "# ERROR plotResults", "$CMGfunc_path/CMGfunc_analyzeResults.pl did not return file $res_filename.func.pdf");	exit }
}
#################################################
#################################################

#=======================================================================================
#=======================================================================================
#===========================	FILETYPE SUBRUTINES	====================================
#=======================================================================================
#=======================================================================================

#################################################
sub checkFileTypeFUNCTIONS {
#################################################
	my $func_filename	= $_[0];
	my $status = "TRUE";
	return $status;
}
#################################################
#################################################

#################################################
sub checkFileTypeRESULTS {
#################################################
	my $status			= "TRUE";
	my $res_filename	= $_[0];
	my $i				= 0;
	#my $fsa_seq_count	= `grep -c '>' $fsa_filename`;
	#my $res_line_count	= `wc $res_filename | awk '{print \$1}'`;
	#my $net_files		= `ls -1 $NET_path | wc | awk '{print \$1}'`;

	# Open the TAB file and test each vector for problems
	open(FP0,$res_filename);
	while( my $line = <FP0> )  {
		next if $line =~ /^#/; 
		my @elements = split("\t", $line);

		unless (scalar(@elements) == 6)	{ 
			printf ("%-40s\t==>\t%-50s%-10s\n", "# ERROR checkFileTypeRESULTS", "RES vectore, $res_filename, line $i does not have the right number of elements");
			print scalar(@elements), "\n";
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
#################################################
#################################################

