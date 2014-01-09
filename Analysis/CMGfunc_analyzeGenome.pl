#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw( sum );
=put
#############################################################################
# Description
#############################################################################

#--------------------------------------------------------------
# Input source
#--------------------------------------------------------------
CMGfunc_testNetworks.py -i file.fsa.tab.vec -n /home/cmgfunc/CMGfunc/NETS/ -m $score_cut

OR

perl CMGfunc.pl -fasta file.fsa -net /home/cmgfunc/CMGfunc/NETS -mse $score_cut

#--------------------------------------------------------------
# Input
#--------------------------------------------------------------

# Comparison started for file:  Escherichia_coli_S88_uid33375_prodigal45d566fb9ab2502b75c554076eaf2117.fsa.tab.vec
Testseq:	13258887 	Network:	proteindata_clan_137.tab.A.mse_0.0269.net.max 	Score:	1.12070115473
Testseq:	13260160 	Network:	proteindata_arch_PF00155_PF00392.tab.C.mse_0.0055.net.max 	Score:	1.10774411783
Testseq:	13260595 	Network:	proteindata_archNOclan_PF00551.tab.C.mse_0.0058.net.max 	Score:	1.03770169883
Testseq:	13261949 	Network:	proteindata_clan_61.tab.B.mse_0.02.net.max 	Score:	1.10672673361
Testseq:	13260249 	Network:	proteindata_archNOclan_PF04055_PF06969.tab.C.mse_0.004.net.max 	Score:	1.00383512348
Testseq:	13258264 	Network:	proteindata_arch_PF00455_PF08220.tab.C.mse_0.0041.net.max 	Score:	1.00314740523
Testseq:	13261714 	Network:	proteindata_clan_48.tab.C.mse_0.0022.net.max 	Score:	0.994838952339

#--------------------------------------------------------------
# Output
#--------------------------------------------------------------
CLAN	425	425	NADH dehydrogenase I, subunit N	1	0.1	425-0008137::0042773::0055114
ARCH	PF00155_PF00392	123:61	Helix-turn-helix clan:PLP dependent aminotransferase superfamily	1	0.1	123:61-0030170::0009058:0003700::0006355::0005622
ARCH	PF03029	NA	NA	2	0.2	PF03029-0000166
ARCH	PF01037_PF01047	123:32	Helix-turn-helix clan:Dimeric alpha/beta barrel superfamily	1	0.1	123:32-0003700::0006355::0005622
ARCH	PF03400	219	Ribonuclease H-like superfamily	2	0.2	219-0003677::0004803::0006313
ARCH	PF00551	NA	NA	5	0.5	PF00551-0016742::0009058
ARCH	PF01514_PF08345	NA	NA	1	0.1	PF01514_PF08345-NA
ARCH	PF02518_PF05231_PF07730	25	His Kinase A (phospho-acceptor) domain	1	0.1	25-0000155::0046983::0000160::0016021:0005524

#--------------------------------------------------------------
# USAGE
#--------------------------------------------------------------
perl CMGfunc_analyzeGenome.pl -res file.fsa.tab.vec.res

cmgfunc@cmgfunc-VirtualBox:~/CMGfunc$ for x in RESTrials/*res
do
perl CMGfunc_analyzeGenome.pl -res $x 
done

=cut

#=======================================================================================
#=======================================================================================
#==============================	ERROR HANDLING	========================================
#=======================================================================================
#=======================================================================================

#.............. Print usage ..............
sub usage {	
	print "#==================================================================================\n";
	printf ("%-40s:\t%-50s\n", "# USAGE", "perl CMGfunc_analyzeGenome.pl -res <CMGfunc res file> -cut <comparison score>");	
	printf ("%-40s:\t%-50s\n", "# USAGE", "-cut is a number between 0 and 1, where 1 would be a perfect match to the fucntional model, default is 0.8");	
	printf ("%-40s:\t%-50s\n", "# EXAMPLE", "perl CMGfunc_analyzeGenome.pl -res file.proteins.fsa.tab.vector.res");
	print "#==================================================================================\n";
	exit( 1 );	}

#.............. Inputs ..............
my ($res, $h, $score_cut) = (undef,undef, 0.8);
my $CMGfunc_path	= "/home/cmgfunc/CMGfunc/indirect";
my $INFO_path		= "/home/cmgfunc/CMGfunc/INFO";

#.............. Get options ..............
&GetOptions ("res:s" =>  \$res, "cut=s" =>  \$score_cut, "h|?" => \$h);
if (defined $h) { &usage }

#.............. Correct input ..............
unless (-d $CMGfunc_path)	{	
	print "#================================================================================================================================\n";	
	printf ("%-40s:\t%-50s\n", "# ERROR checkPATHS", "nothing called $CMGfunc_path or is not a directory"); 
	printf ("%-40s:\t%-50s\n", "# ERROR checkPATHS", "CMGfunc_directory : $CMGfunc_path");		
	&usage;	exit; }

unless ($res)	{	
	print "#================================================================================================================================\n";	
	printf ("%-40s:\t%-50s\n", "# ERROR", "CMGfunc RES file not defined, example: perl CMGfunc_analyzeGenome.pl -res file.fsa.tab.vector.res");	
	&usage;	exit; }

unless ($score_cut >=0 and $score_cut <= 1) {
	print "#================================================================================================================================\n";	
	printf ("%-40s:\t%-50s\n", "# ERROR", "Cutoff score is not a number between 0 and 1: $score_cut");	
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

if (-e $res.".table") { printf ("%-40s\t==>\t%-50s\n", "# FILE", "already exists : $res.table") }
else {	&analyzeResults($res)	}	#.............. Analyze results ..............

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
	my (%GOHash, %freqHash);

	printf ("%-40s\t==>\t%-50s\n", "# RUNNING analyzeResults", "will create file: $res_filename.table ");

	open(FHres,$res_filename);
	while( my $line = <FHres> )  {
        chomp $line;

		next if $line =~ m/^#/ or $line =~ m/group_concat/;
		my @elements = split("\t", $line);

		my ($seq, $net, $score) = ($elements[1], $elements[3], $elements[5]);
		
		if ($score >= $score_cut) {
			if		($net =~ m/\_clan\_/) 		{	&clanNetworks		($net, $score, \%GOHash, \%freqHash)	}
			elsif	($net =~ m/\_arch\_/)		{	&archClanNetworks	($net, $score, \%GOHash, \%freqHash)	}
			elsif	($net =~ m/\_archNOclan\_/)	{	&archNoClanNetworks	($net, $score, \%GOHash, \%freqHash)	}
			elsif	($net =~ m/local/)			{	&localNetworks		($net, $score, \%freqHash)				}
		}
	} 
	close FHres;# WHILE RESLINE END  ---------------------------------------------

	my $sum = 0;
	foreach my $k (keys %freqHash) { $sum += $freqHash{$k}}

	open(FHtable, ">", $res_filename.".table");	
	print FHtable "#================================================================================================\n";
	print FHtable "# Total number of proteins in set :\t",$res_seq_count,"\n";	
	print FHtable "# Cutoff for acceptable match between protein and function :\t",$score_cut,"\n";	
	print FHtable "# Proteins with matches to functional models :\t",$sum,"\n";	
	print FHtable "# Proteins without matches to functional models :\t",$res_seq_count-$sum,"\n";	
	print FHtable "# Type\tmodel nr.\tclans\tdescription\tfreq\tpercentage\tgo terms\n";
	print FHtable "#================================================================================================\n";
	
	foreach my $n ( keys %freqHash ) {
		if ($n =~ m/local/) {
			my $netnr = (split("\\.", $n))[0];	
			my $desc = `head -n 1 /home/cmgfunc/CMGfunc/NETS/LOCAL/$net*txt | gut -f 2`;
			my $name = `head -n 1 /home/cmgfunc/CMGfunc/NETS/LOCAL/$net*txt | gut -f 1`;	

			print FHtable "LOCAL\t$n\t$n\t$name\t$freqHash{$n}\tNA\t$desc\tNA\tNA\n";
		}
		if ($n =~ m/PF/) {
			my ($pfam, $clan_nr, $clan_name, $pfam_desc, $go_terms);	
			my (%clan_nrs, %clan_names, %pfam_descs, %go_terms_all);
			my @pfams = split("_", $n);
			my (@go_mf, @go_cc, @go_bp);

			for $pfam (@pfams) {
				$clan_nr = `grep -P -m 1 "^$pfam\t" /home/cmgfunc/CMGfunc/INFO/info_pfam_clannr_clandesc_gos_godescs.txt | cut -f 3`; chomp $clan_nr;
				$clan_nrs{$clan_nr} = 1 unless $clan_nr eq "NA";

				$clan_name = `grep -P -m 1 "^$pfam\t" /home/cmgfunc/CMGfunc/INFO/info_pfam_clannr_clandesc_gos_godescs.txt | cut -f 4`; chomp $clan_name;
				$clan_names{$clan_name} = 1 unless $clan_name eq "NA";

				$pfam_desc = `grep -P -m 1 "^$pfam\t" /home/cmgfunc/CMGfunc/INFO/info_pfam_clannr_clandesc_gos_godescs.txt | cut -f 2`; chomp $pfam_desc;
				$pfam_descs{$pfam_desc} = 1 unless $pfam_desc eq "NA";

				$go_terms = `grep -P -m 1 "^$pfam\t" /home/cmgfunc/CMGfunc/INFO/info_pfam_clannr_clandesc_gos_godescs.txt | cut -f 6`; chomp $go_terms;
				my @go_tmp = split("::", $go_terms);
				for my $g (@go_tmp) {
					next if $g eq "NA";
					my ($go_type, $type_go) = ('NA', 'NA');
					$go_type = `grep -P $g /home/cmgfunc/CMGfunc/INFO/info_goid_desc_type_go.txt | cut -f 3`; chomp $go_type;
					if		($go_type =~ m/molecular/)	{	$type_go = $g."MF"; push(@go_mf,$type_go) } 
					elsif	($go_type =~ m/cellular/)	{	$type_go = $g."CC"; push(@go_cc,$type_go) } 
					elsif	($go_type =~ m/biological/) {	$type_go = $g."BP"; push(@go_bp,$type_go) }
				}		

				$go_terms_all{join(":", @go_mf)} = 1 unless $go_terms eq "NA";
				$go_terms_all{join(":", @go_cc)} = 1 unless $go_terms eq "NA";
				$go_terms_all{join(":", @go_bp)} = 1 unless $go_terms eq "NA";
			}
			my ($clan_nrs_line, $clan_names_line, $go_terms_line, $pfam_descs_line, $name_line) = ("NA","NA","NA","NA","NA");
			$clan_nrs_line = join(":", (keys(%clan_nrs))) 		unless scalar(keys(%clan_nrs)) < 1;
			$clan_names_line = join(":", (keys(%clan_names))) 	unless scalar(keys(%clan_names)) < 1;
			$go_terms_line = join(":", (keys(%go_terms_all))) 	unless scalar(keys(%go_terms_all)) < 1;
			$pfam_descs_line = join(":", (keys(%pfam_descs))) 	unless scalar(keys(%pfam_descs)) < 1;
			
			#$name_line =  $clan_nrs_line."-".$go_terms_line if scalar(keys(%clan_nrs)) > 0 ;
			#$name_line =  $n."-".$go_terms_line if scalar(keys(%clan_nrs)) < 1 ;
			my ($go_terms_mf, $go_terms_cc,$go_terms_bp) = ("NA","NA","NA");
			if (scalar(@go_mf) > 0) { $go_terms_mf = join(":", @go_mf); $go_terms_mf =~ s/::$//g }
			if (scalar(@go_cc) > 0) { $go_terms_cc = join(":", @go_cc); $go_terms_cc =~ s/::$//g }
			if (scalar(@go_bp) > 0) { $go_terms_bp = join(":", @go_bp); $go_terms_bp =~ s/::$//g }
			my $name_line_mf =  $n."-".$go_terms_mf;
			my $name_line_cc =  $n."-".$go_terms_cc;
			my $name_line_bp =  $n."-".$go_terms_bp;

			my $percentage = ($freqHash{$n}/$res_seq_count)*100;

			print FHtable "ARCH\t$n\t$clan_nrs_line\t$clan_names_line\t$freqHash{$n}\t$percentage\t$name_line_mf\t$name_line_cc\t$name_line_bp\n";
			#print FHtable "ARCH\t$n\t$clan_nrs_line\t$clan_names_line\t$freqHash{$n}\t$percentage\t$go_terms_line\n" if scalar(keys(%clan_nrs)) > 0;
			#print FHtable "ARCH\t$n\tNA\t$pfam_descs_line\t$freqHash{$n}\t$percentage\t$go_terms_line\n" if scalar(keys(%clan_nrs)) < 1;
		}
		else {
			my (@go_mf, @go_cc, @go_bp, @go_tmp2);
			my %go_terms_all;
			my ($go_type, $type_go) = ('NA', 'NA');

			my $clan_name = `cut -f 3,4 /home/cmgfunc/CMGfunc/INFO/info_pfam_clannr_clandesc_gos_godescs.txt | grep -P "^$n\t" | cut -f 2 | sort -u`; chomp $clan_name;
			my $go_terms = `cut -f 3,6 /home/cmgfunc/CMGfunc/INFO/info_pfam_clannr_clandesc_gos_godescs.txt | grep -P "^$n\t" | cut -f 2 | perl -pe 's/::/\n/g' | sort -u | 				grep -v NA | perl -pe 's/\n/::/g'`;
			my @go_tmp = split("::", $go_terms);
			
			for my $g (@go_tmp) {
				next if $g eq "NA";
				my ($go_type, $type_go) = ('NA', 'NA');
				$go_type = `grep -P $g /home/cmgfunc/CMGfunc/INFO/info_goid_desc_type_go.txt | cut -f 3`; chomp $go_type;
				if		($go_type =~ m/molecular/)	{	$type_go = $g."MF"; push(@go_mf,$type_go) } 
				elsif	($go_type =~ m/cellular/)	{	$type_go = $g."CC"; push(@go_cc,$type_go) } 
				elsif	($go_type =~ m/biological/) {	$type_go = $g."BP"; push(@go_bp,$type_go) }
			}		

			$go_terms_all{join(":", @go_mf)} = 1 unless $go_terms eq "NA";
			$go_terms_all{join(":", @go_cc)} = 1 unless $go_terms eq "NA";
			$go_terms_all{join(":", @go_bp)} = 1 unless $go_terms eq "NA";

			my ($go_terms_mf, $go_terms_cc,$go_terms_bp) = ("NA","NA","NA");
			if (scalar(@go_mf) > 0) { $go_terms_mf = join(":", @go_mf); $go_terms_mf =~ s/::$//g }
			if (scalar(@go_cc) > 0) { $go_terms_cc = join(":", @go_cc); $go_terms_cc =~ s/::$//g }
			if (scalar(@go_bp) > 0) { $go_terms_bp = join(":", @go_bp); $go_terms_bp =~ s/::$//g }
			my $name_line_mf =  $n."-".$go_terms_mf;
			my $name_line_cc =  $n."-".$go_terms_cc;
			my $name_line_bp =  $n."-".$go_terms_bp;


			my $percentage = ($freqHash{$n}/$res_seq_count)*100;


			print FHtable "CLAN\t$n\t$n\t$clan_name\t$freqHash{$n}\t$percentage\t$name_line_mf\t$name_line_cc\t$name_line_bp\n";
		}
	}
	close FHtable;
}	# SUB RUTINE END  ---------------------------------------------

#=======================================================================================
#=======================================================================================
#==========================	CALCULATION SUBRUTINES	====================================
#=======================================================================================
#=======================================================================================

#################################################
sub localNetworks {
#################################################
	my ($net, $score, $freqHash) = @_;
	$$freqHash{$net} ++ ;
}

#################################################
sub clanNetworks {
#################################################
	my ($net, $score, $GOHash, $freqHash) = @_;

	# HASH of clan networks	---------------------------------------------
	my $clannr = (split("_",(split("\\.", $net))[0]))[-1];	

	$$freqHash{$clannr} ++;
}

#################################################
#################################################

#################################################
sub archClanNetworks {
#################################################
	my ($net, $score, $GOHash, $freqHash) = @_;

	# HASH of Archclan networks	-------------------------------------
	my $arch = (split("\\.", $net))[0]; $arch =~ s/proteindata_arch_//g;
	my @pfams = (split("_", $arch));
	my %go_per_arch;

	# Seperate high and bad scoring results --------------------------------
	$$freqHash{$arch} ++ ;
}
#################################################
#################################################

#################################################
sub archNoClanNetworks {
#################################################
	my ($net, $score, $GOHash, $freqHash) = @_;

	# HASH of ArchNoclan networks	---------------------------------
	my $arch = (split("\\.", $net))[0]; $arch =~ s/proteindata_archNOclan_//g;
	my @pfams = (split("_", $arch));
	my %go_per_arch;

	$$freqHash{$arch} ++ ;
}
#################################################
#################################################

#=======================================================================================
#=======================================================================================
#===========================	FILETYPE SUBRUTINES	====================================
#=======================================================================================
#=======================================================================================

#################################################
#################################################

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
#################################################
#################################################

