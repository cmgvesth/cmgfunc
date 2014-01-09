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
perl CMGfunc_calcFeatures.pl -fasta genomeA.proteins.fsa

#--------------------------------------------------------------
# Input
#--------------------------------------------------------------
head trial.fsa.tab
protein_id	protein_id	fasta	absorbtion_max	absorbtion_min	aliphatic_index	amino_acid_percentages_A	amino_acid_percentages_C	amino_acid_percentages_D	amino_acid_percentages_E	amino_acid_percentages_F	amino_acid_percentages_G	amino_acid_percentages_H	amino_acid_percentages_I	amino_acid_percentages_K	amino_acid_percentages_L	amino_acid_percentages_M	amino_acid_percentages_N	amino_acid_percentages_P	amino_acid_percentages_Q	amino_acid_percentages_R	amino_acid_percentages_S	amino_acid_percentages_T	amino_acid_percentages_V	amino_acid_percentages_W	amino_acid_percentages_Y	aromaticity	extinction_coeff_maxextinction_coeff_min	flexibility_average	gravy_coeff	in_vivo_halflife_ecoli	in_vivo_halflife_human	in_vivo_halflife_yeast	instability_index	isoelectric_point	molecular_weight	number_aromatic	number_negative	number_nonpolar	number_positive	number_uncharged	secondary_structure_helix	secondary_structure_sheet	secondary_structure_turn	sequence_length	total_atoms	total_molecules	weight_aromatic	weight_negative	weight_nonpolar	weight_positive	weight_uncharged	Cytoplasmic_Score_n	CytoplasmicMembrane_Score_n	Periplasmic_Score_n	OuterMembrane_Score_n	Extracellular_Score_n	Final_Score_n	Cytoplasmic_Score_p	CytoplasmicMembrane_Score_p	Cellwall_Score_p	Extracellular_Score_p	Final_Score_p	netphos_total	netphos_serine	netphos_threonine	Cmax_n	Ymax_n	Smax_n	Smean_n	D_n	signalp_n	Cmax_p	Ymax_p	Smax_p	Smean_p	D_p	signalp_p	ExpAA	First60	PredHel	seq_low_total	seq_low_avg	seq_low_coverage	seq_high_total	seq_high_avg	seq_high_coverage
9832486_PF00004_23	9832486_PF00004_23	trial.fsa	0.002	0.002	84.574	0.134	0.005	0.061	0.093	0.033	0.067	0.018	0.025	0.038	0.111	0.031	0.031	0.043	0.021	0.074	0.044	0.067	0.062	0.011	0.030	0.074	107.389	107.082	1.001	-0.285	600	6000	1200	45.160	5.008	67286.210	0.074	0.154	0.523	0.130	0.193	0.272	0.370	0.185	610	18.428	610	7994.830	13311.110	36648.150	12908.130	15393	8.96	0.51	0.26	0.01	0.26	8.96	7.501.15	0.62	0.73	7.50	1079	390	625	0.115	0.175	0.414	0.275	0.222	0	0.142	0.221	0.501	0.290	0.248	0	0.60	0.59	0	3	16.6666666666667	2.73224043715847	4	138.25	22.6639344262295


#--------------------------------------------------------------
# Output
#--------------------------------------------------------------
head trial.fsa.tab.vec
9272486_PF00005_PF00664_23	 trial.fsa	 9272486_PF00005_PF00664_23	 0	 0	 0	 1	 1	 0	 0	 1	 1	 0	 0	 0.232817869415808	 0.232817869415808	 0.00687285223367698	 0.0189003436426117	 0.0189003436426117	 0.00515463917525773	 0.122	 0.118	 0.143	 0.176	 0	 0	 0.243	 0.373	 0.134	 0.179	 0.148	 0.174	 0.024404	 0.373666666666667	 0.0166666666666667	 0.001	 0.001	 0.2571725	 0.084	 0.003	 0.052	 0.052	 0.040	 0.064	 0.015	 0.074	 0.043	 0.108	 0.046	 0.036	 0.021	 0.045	 0.062	 0.086	 0.058	 0.081	 0.009	 0.022	 0.000120274914089347	 0.0448594444444444	 0.0447405555555556	 0.903636363636364	 0.1	 0.170916	 0.553714285714286	 0.0136898461303436	 0.070	 0.103	 0.529	 0.120	 0.247	 0.333	 0.290	 0.206	 0.01455	 0.0321769759450172	 0.111218861852528	 0.130296425957675	 0.573875248096211	 0.175482871248489	 0.282611458472499	 
7242486_PF00005_PF08352_23	 trial.fsa	 7242486_PF00005_PF08352_23	 0.008	 0.211	 0.105	 0.788	 0.878	 0	 0.009	 0.788	 0.878	 0	 0	 0.475460122699386	 0.475460122699387	 0.00613496932515337	 0.0398773006134969	 0.0398773006134969	 0.00306748466257669	 0.152	 0.153	 0.114	 0.183	 0	 0	 0.175	 0.240	 0.1

#--------------------------------------------------------------
# USAGE
#--------------------------------------------------------------
perl CMGfunc_normFeatures.pl -tab somefile.tab

Creates file.fsa.tab.vec
=cut

my $filename;
&GetOptions ("tab:s" =>  \$filename);

#.............. Input files and paths exists ..............
unless (defined $filename)		{	
	printf ("%-30s:\t%-50s\n", "# ERROR", "Tab file not defined");	
	printf ("%-30s:\t%-50s\n", "# USAGE", "perl CMGfunc_normFeatures.pl -tab somefile.tab");	
	exit; 
}

#------------------------------------------
#-------- Vectores from file --------------
#------------------------------------------
my (@AoA, @filenames, @tmp, %HoH);

my $linecount=0;
open(FILE,$filename);

while( my $line = <FILE> )  {
    chomp $line;
	if ($line =~ m/number/) {	@filenames = split('\t', $line)	}
    else 					{
		@tmp = split('\t', $line); $linecount++;
		for my $i ( 0 .. $#filenames )	{	
			$filenames[$i] =~ s/\s//;
			$HoH{$linecount}{$filenames[$i]} = $tmp[$i];	
		}
	}
} 
close FILE;

#print Dumper(@filenames);
#foreach my $f ( sort keys %HoH ) {	for my $r ( sort keys %{ $HoH{$f} } ) {	print "$f: $r=$HoH{$f}{$r} \n"	}	}
#for my $i ( 0 .. $#filenames ) {	print "$filenames[$i] :: $i\n"	}

#------------------------------------------
#-------- Feature names per column --------
#------------------------------------------
my $name_line = "protein_id uarch auto_clan nCellwall_Score_p nCytoplasmic_Score_n nCytoplasmic_Score_p nCytoplasmicMembrane_Score_n nCytoplasmicMembrane_Score_p nExtracellular_Score_n nExtracellular_Score_p nFinal_Score_n nFinal_Score_p nOuterMembrane_Score_n nPeriplasmic_Score_n nseq_high_avg nseq_high_coverage nseq_high_total nseq_low_avg nseq_low_coverage nseq_low_total Cmax_n Cmax_p D_n D_p signalp_n signalp_p Smax_n Smax_p Smean_n Smean_p Ymax_n Ymax_p absorbtion_max absorbtion_min naliphatic_index amino_acid_percentages_A amino_acid_percentages_C amino_acid_percentages_D amino_acid_percentages_E amino_acid_percentages_F amino_acid_percentages_G amino_acid_percentages_H amino_acid_percentages_I amino_acid_percentages_K amino_acid_percentages_L amino_acid_percentages_M amino_acid_percentages_N amino_acid_percentages_P amino_acid_percentages_Q amino_acid_percentages_R amino_acid_percentages_S amino_acid_percentages_T amino_acid_percentages_V amino_acid_percentages_W amino_acid_percentages_Y naromaticity nextinction_coeff_max nextinction_coeff_min nflexibility_average nin_vivo_halflife_ecoli ninstability_index nisoelectric_point nmolecular_weight number_aromatic number_negative number_nonpolar number_positive number_uncharged secondary_structure_helix secondary_structure_sheet secondary_structure_turn nsequence_length ntotal_atoms nweight_aromatic nweight_negative nweight_nonpolar nweight_positive nweight_uncharged";

my @name_array = split('\s', $name_line);
my %nameHash;
foreach my $item (@name_array) { $nameHash{$item} = 1 }

my %normHoH;

foreach my $f ( sort keys %HoH ) {	
	foreach my $r ( sort keys %{ $HoH{$f} } ) {	
		
		foreach my $i ( keys %nameHash ) {
			if ($i =~ m/$r/) {
				#print "$f :: $r :: $HoH{$f}{$r} :: $i\n";
				if 		($r =~ m/aliphatic_index/)			{	$normHoH{$f}{$i} = $HoH{$f}{$r}/400	}
				elsif	($r =~ m/aromaticity/)				{	$normHoH{$f}{$i} = $HoH{$f}{$r}/$HoH{$f}{"sequence_length"}	}
				elsif	($r =~ m/extinction_coeff_max/)		{	$normHoH{$f}{$i} = $HoH{$f}{$r}/1800	}
				elsif	($r =~ m/extinction_coeff_min/)		{	$normHoH{$f}{$i} = $HoH{$f}{$r}/1800	}
				elsif	($r =~ m/in_vivo_halflife_ecoli/)	{	$normHoH{$f}{$i} = $HoH{$f}{$r}/6000	}
				elsif	($r =~ m/instability_index/)		{	$normHoH{$f}{$i} = $HoH{$f}{$r}/250	}
				elsif	($r =~ m/total_atoms/)				{	$normHoH{$f}{$i} = $HoH{$f}{$r}/$HoH{$f}{"sequence_length"}	}
				elsif	($r =~ m/isoelectric_point/)		{	$normHoH{$f}{$i} = $HoH{$f}{$r}/14	}
				elsif	($r =~ m/molecular_weight/)			{	$normHoH{$f}{$i} = $HoH{$f}{$r}/4713080	}
				elsif	($r =~ m/sequence_length/)			{	$normHoH{$f}{$i} = $HoH{$f}{$r}/40000	}
				elsif	($r =~ m/weight_aromatic/)			{	$normHoH{$f}{$i} = $HoH{$f}{$r}/$HoH{$f}{"molecular_weight"}	}
				elsif	($r =~ m/weight_negative/)			{	$normHoH{$f}{$i} = $HoH{$f}{$r}/$HoH{$f}{"molecular_weight"}	}
				elsif	($r =~ m/weight_nonpolar/)			{	$normHoH{$f}{$i} = $HoH{$f}{$r}/$HoH{$f}{"molecular_weight"}	}
				elsif	($r =~ m/weight_positive/)			{	$normHoH{$f}{$i} = $HoH{$f}{$r}/$HoH{$f}{"molecular_weight"}	}
				elsif	($r =~ m/weight_uncharged/)			{	$normHoH{$f}{$i} = $HoH{$f}{$r}/$HoH{$f}{"molecular_weight"}	}
				#elsif	($r =~ m/ExpAA/)					{	$normHoH{$f}{$i} = $HoH{$f}{$r}/5000	}
				#elsif	($r =~ m/First60/)					{	$normHoH{$f}{$i} = $HoH{$f}{$r}/60	}
				#elsif	($r =~ m/PredHel/)					{	$normHoH{$f}{$i} = $HoH{$f}{$r}/300	}
				elsif	($r =~ m/flexibility_average/)		{	$normHoH{$f}{$i} = $HoH{$f}{$r}/1.1	}
				elsif	($r =~ m/Score_p/)					{	$normHoH{$f}{$i} = $HoH{$f}{$r}/10	}
				elsif	($r =~ m/Score_n/)					{	$normHoH{$f}{$i} = $HoH{$f}{$r}/10	}
				elsif	($r =~ m/seq_high_avg/)				{	$normHoH{$f}{$i} = $HoH{$f}{$r}/$HoH{$f}{"sequence_length"}	}
				elsif	($r =~ m/seq_high_total/)			{	$normHoH{$f}{$i} = $HoH{$f}{$r}/$HoH{$f}{"sequence_length"}	}
				elsif	($r =~ m/seq_high_coverage/)		{	$normHoH{$f}{$i} = $HoH{$f}{$r}/100	}
				elsif	($r =~ m/seq_low_avg/)				{	$normHoH{$f}{$i} = $HoH{$f}{$r}/$HoH{$f}{"sequence_length"}	}
				elsif	($r =~ m/seq_low_total/)			{	$normHoH{$f}{$i} = $HoH{$f}{$r}/$HoH{$f}{"sequence_length"}	}
				elsif	($r =~ m/seq_low_coverage/)			{	$normHoH{$f}{$i} = $HoH{$f}{$r}/100	}
				else 										{	$normHoH{$f}{$i} = $HoH{$f}{$r}	}	
			}
				
			elsif	($i =~ m/protein_id/)				{	$normHoH{$f}{$i} = $HoH{$f}{"protein_id"}	}
			elsif	($i =~ m/uarch/)					{	$normHoH{$f}{$i} = $HoH{$f}{"fasta"}	}
			elsif	($i =~ m/auto_clan/)				{	$normHoH{$f}{$i} = $HoH{$f}{"protein_nr"}	}
		}
	}	
}
open (OUTFILE, ">", $filename.".vec") or die printf ("%-40s\t==>\t%-50s\n", "# ERROR", "cannot create out vector file, $filename.vec");
#print Dumper(%normHoH);

foreach my $m ( keys %normHoH ) {
	open (OUTFILE, ">>", $filename.".vec") or die printf ("%-40s\t==>\t%-50s\n", "# ERROR", "cannot add to output vector file, $filename.vec");
	foreach my $i ( @name_array ) {	
		#print "$i :: $m :: $normHoH{$m}{$i}\n";
		print OUTFILE "$normHoH{$m}{$i}\t "
		#print "$m: $i=$normHoH{$m}{$i} \n"
	}
	print OUTFILE "\n";
	`sed -i 's/N/0/' $filename.vec`;
	close OUTFILE ;
}


