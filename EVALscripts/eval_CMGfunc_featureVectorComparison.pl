#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Data::Dumper;
use Bio::SeqIO;

=put

#############################################################################
# Feature vector comparison
#############################################################################

#--------------------------------------------------------------
# Input source
#--------------------------------------------------------------
# Split data based on clans (takes 20 min):
awk  '{print >  "proteindata_clan_"$3".tab"}' data_feat_all_clan_unique.txt

# Split large clans on architectures:
for x in proteindata_clan_*tab; do LINES=$(grep -c -v "uarch" $x | cut -f1) ; if [ $LINES -gt 10000 ]; then awk  '{print >  "proteindata_arch_"$2".tab"}' $x; mv $x $x.large ; fi; done

#--------------------------------------------------------------
# Input
#--------------------------------------------------------------
head proteindata_clan_101.tab 
10002112	PF00588_PF08032	101	0.25	0.926	0.25	0.024	0.25	0.001	0.25	0.926	0.25	0.001	0.048	0.996	0.996	0.004	0	0	0	0.103	0.109	0.117	0.128	00	0.152	0.204	0.123	0.136	0.111	0.122	0.000416	0	0	0.004	0.004	0.2331	0.104	0.012	0.06	0.064	0.016	0.076	0.032	0.052	0.044	0.092	0.036	0.064	0.036	0.02	0.064	0.052	0.048	0.092	0.012	0.024	0.000208	0.05695	0.0565333	0.907273	0.1	0.101632	0.4335	0.00581429	0.052	0.124	0.528	0.14	0.208	0.288	0.296	0.228	0.00625	0.073408	0.0861429	0.158762	0.553789	0.20569	0.245498
10004027	PF00588_PF08032	101	0.062	0.926	0.75	0.024	0.115	0.001	0.073	0.926	0.75	0.001	0.048	0.996	0.996	0.004	0	0	0	0.102	0.108	0.115	0.128	00	0.15	0.204	0.121	0.137	0.11	0.122	0.000402	0	0	0.004	0.004	0.2331	0.1	0.016	0.056	0.068	0.016	0.072	0.032	0.056	0.044	0.092	0.036	0.06	0.036	0.02	0.064	0.06	0.048	0.088	0.012	0.024	0.000208	0.0570889	0.0565333	0.907273	0.1	0.113872	0.433643	0.00582768	0.052	0.124	0.524	0.14	0.212	0.288	0.296	0.228	0.00625	0.07352	0.085945	0.158908	0.551461	0.205218	0.247776

#--------------------------------------------------------------
# Output
#--------------------------------------------------------------
perl /home/tammi/MySQL/CMGfunc/eval_CMGfunc_featureVectorComparison.pl -vec CLANS/proteindata_clan_101.tab
nCytoplasmic_Score_n          	AVE: 0.3926	STD: 0.3164
nCytoplasmic_Score_p          	AVE: 0.4648	STD: 0.3022
nCytoplasmicMembrane_Score_n  	AVE: 0.1567	STD: 0.0712
nCytoplasmicMembrane_Score_p  	AVE: 0.2701	STD: 0.2664
nExtracellular_Score_n        	AVE: 0.1503	STD: 0.0813
nExtracellular_Score_p        	AVE: 0.1339	STD: 0.1032
nFinal_Score_n                	AVE: 0.3930	STD: 0.3161
nFinal_Score_p                	AVE: 0.5774	STD: 0.2920
nOuterMembrane_Score_n        	AVE: 0.1455	STD: 0.0889
nPeriplasmic_Score_n          	AVE: 0.1535	STD: 0.0761
nseq_high_avg                 	AVE: 0.8760	STD: 0.2364
nseq_high_coverage            	AVE: 0.8759	STD: 0.2364
nseq_high_total               	AVE: 0.0119	STD: 0.0044
nseq_low_avg                  	AVE: 0.0228	STD: 0.0493
nseq_low_coverage             	AVE: 0.0208	STD: 0.0466
nseq_low_total                	AVE: 0.0021	STD: 0.0045
Cmax_n                        	AVE: 0.1210	STD: 0.0225
Cmax_p                        	AVE: 0.1380	STD: 0.0447
D_n                           	AVE: 0.1383	STD: 0.0457
.....

#--------------------------------------------------------------
# USAGE
#--------------------------------------------------------------
for x in proteindata_*tab
do
LINES=$(grep -c -v "uarch" $x | cut -f1) 
if [ $LINES -gt 10000 ]; then
	echo "# To many seqs in $x, nr. lines $LINES" 
	continue;
fi
perl /home/tammi/MySQL/CMGfunc/eval_CMGfunc_featureVectorComparison.pl -vec $x > $x.veccomp
done

# Check no values bend the 0-1 range
grep -v "STD: 0\." *comp
grep -v "AVE: 0\." *comp


=cut

my $vectorfilename;
&GetOptions ("vec:s" =>  \$vectorfilename);

#.............. Input files and paths exists ..............
unless (defined $vectorfilename)		{	
	printf ("%-30s:\t%-50s\n", "# ERROR", "Vector file not defined");	
	printf ("%-30s:\t%-50s\n", "# USAGE", "perl SCRIPTS/compareVector.pl -vec somefile.tab");	
	exit; }

#------------------------------------------
#-------- Vectores from file --------------
#------------------------------------------

my @AoA;

open(FILE,$vectorfilename);
while( my $line = <FILE> )  {
    chomp $line;
	next if $line =~ m/\_/;
    my @tmp = split('\t', $line);
    push @AoA, [ @tmp[0 .. $#tmp] ];
} 
#------------------------------------------
#-------- Hash of array of vectores -------
#------------------------------------------

my %HoA;

for my $i ( 0 .. $#AoA ) {
        for my $j ( 0 .. $#{ $AoA[$i] }  ) {
			push(@{$HoA{$j}}, $AoA[$i][$j]);
    }
}

#------------------------------------------
#-------- Feature names per column --------
#------------------------------------------

my $raw_line = "protein_id	uarch	auto_clan	nCellwall_Score_p	nCytoplasmic_Score_n	nCytoplasmic_Score_p	nCytoplasmicMembrane_Score_n	nCytoplasmicMembrane_Score_p	nExtracellular_Score_n	nExtracellular_Score_p	nFinal_Score_n	nFinal_Score_p	nOuterMembrane_Score_n	nPeriplasmic_Score_n	nseq_high_avg	nseq_high_coverage	nseq_high_total	nseq_low_avg	nseq_low_coverage	nseq_low_total	Cmax_n	Cmax_p	D_n	D_p	signalp_n	signalp_p	Smax_n	Smax_p	Smean_n	Smean_p	Ymax_n	Ymax_p	ExpAA	First60	PredHel	absorbtion_max	absorbtion_min	naliphatic_index	amino_acid_percentages_A	amino_acid_percentages_C	amino_acid_percentages_D	amino_acid_percentages_E	amino_acid_percentages_F	amino_acid_percentages_G	amino_acid_percentages_H	amino_acid_percentages_I	amino_acid_percentages_K	amino_acid_percentages_L	amino_acid_percentages_M	amino_acid_percentages_N	amino_acid_percentages_P	amino_acid_percentages_Q	amino_acid_percentages_R	amino_acid_percentages_S	amino_acid_percentages_T	amino_acid_percentages_V	amino_acid_percentages_W	amino_acid_percentages_Y	naromaticity	nextinction_coeff_max	nextinction_coeff_min	flexibility_average	nin_vivo_halflife_ecoli	ninstability_index	nisoelectric_point	nmolecular_weight	number_aromatic	number_negative	number_nonpolar	number_positive	number_uncharged	secondary_structure_helix	secondary_structure_sheet	secondary_structure_turn	nsequence_length	ntotal_atoms	nweight_aromatic	nweight_negative	nweight_nonpolar	nweight_positive	nweight_uncharged";

my @raw_array = split('\t', $raw_line);
my $count = 0;
my %NameHash;

foreach my $name (@raw_array) {
	$NameHash{$count} = $name;
	$count++;
}

#------------------------------------------
#-------- AVG and STD per feature ---------
#------------------------------------------
my ($ave, $std) = ("", "");

foreach my $element ( sort {$a <=> $b} (keys %HoA) ) {
	#print @{ $HoA{$element} }, "\n";
	$ave = &average(\@{ $HoA{$element} }) unless $element < 3;
	$std = &stdev(\@{ $HoA{$element} }) unless $element < 3;
	printf ("%-30s\t%-5s%.4f\t%-5s%.4f\n", "$NameHash{$element}", "AVE:", $ave, "STD:", $std) if $element > 3;	
	#printf ("%-30s\t%-5s%.4f\t%-5s%\s\n", "$NameHash{$element}", " ", $element, "STD:", $element) if $element <= 3;	
	#printf "$NameHash{$element}\t\t\tAVE:\t%.4f\t\tSTD:\t%.4f\n", $ave, $std;
}

#------------------------------------------
#------------------ SUBS ------------------
#------------------------------------------
sub average	{
        my($data) = @_;
        if (not @$data) {	die("Empty array\n")	}
        my $total = 0;
        foreach (@$data) {	$total += $_	}
        my $average = $total / @$data;
        return $average;
}
sub stdev	{
        my($data) = @_;
		if(@$data == 1){	return 0	}
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {	$sqtotal += ($average-$_) ** 2	}
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

