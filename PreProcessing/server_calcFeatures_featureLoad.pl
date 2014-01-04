#!/tools/bin/perl

use DBI;
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Digest::MD5 qw(md5 md5_hex);

=put
#------------------------------------------------------------------------
# Doc
# http://www.cbs.dtu.dk/services/ProtFun/protfun_add.php

CREATE TABLE feat (
protein_id int(11) DEFAULT NULL, absorbtion_max float, absorbtion_min float, 
aliphatic_index float, amino_acid_percentages_A float, amino_acid_percentages_C float, 
amino_acid_percentages_D float, amino_acid_percentages_E float, amino_acid_percentages_F float, 
amino_acid_percentages_G float, amino_acid_percentages_H float, amino_acid_percentages_I float, 
amino_acid_percentages_K float, amino_acid_percentages_L float, amino_acid_percentages_M float, 
amino_acid_percentages_N float, amino_acid_percentages_P float, amino_acid_percentages_Q float, 
amino_acid_percentages_R float, amino_acid_percentages_S float, amino_acid_percentages_T float, 
amino_acid_percentages_V float, amino_acid_percentages_W float, amino_acid_percentages_Y float, 
aromaticity float, extinction_coeff_max float, extinction_coeff_min float, flexibility_average float, 
gravy_coeff float, in_vivo_halflife_ecoli float, in_vivo_halflife_human float, in_vivo_halflife_yeast float, 
instability_index float, isoelectric_point float, molecular_weight float, number_aromatic float, 
number_negative float, number_nonpolar float, number_positive float, number_uncharged float, 
secondary_structure_helix float, secondary_structure_sheet float, secondary_structure_turn float, 
sequence_length float, total_atoms float, total_molecules float, weight_aromatic float, weight_negative float, 
weight_nonpolar float, weight_positive float, weight_uncharged float, Cytoplasmic_Score_n float, 
CytoplasmicMembrane_Score_n float, Periplasmic_Score_n float, OuterMembrane_Score_n float, Extracellular_Score_n float, 
Final_Score_n float, Cytoplasmic_Score_p float, CytoplasmicMembrane_Score_p float, Cellwall_Score_p float, 
Extracellular_Score_p float, Final_Score_p float, netphos_total float, netphos_serine float, netphos_threonine float, 
Cmax_n float, Ymax_n float, Smax_n float, Smean_n float, D_n float, signalp_n float, Cmax_p float, Ymax_p float, 
Smax_p float, Smean_p float, D_p float, signalp_p float, ExpAA float, First60 float, PredHel float, seq_low_total float, 
seq_low_avg float, seq_low_coverage float, seq_high_total float, seq_high_avg float, seq_high_coverage float, 
UNIQUE KEY uk_protein_id (protein_id), CONSTRAINT fk_protein_id_feat FOREIGN KEY (protein_id) REFERENCES protein (protein_id) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
=cut
#------------------------------------------------------------------------

#=============================================================================================
#=============================== Set defaults ===================================
#=============================================================================================

my $dbuser = '<insert user name>';
my $dbpass = '<insert password>';
my $dbstring = '<insert dbi path>';

my $analysis_directory = "/home/cmgfunc/CMGfunc/ProtFun";


#=============================================================================================
#=============================== Set database ===================================
#=============================================================================================

#TODO: This needs to be changed to your actual settings on the server
my $dbh = DBI-> connect($dbstring,$dbuser,$dbpass, 
	{ RaiseError => 0, AutoCommit => 0, PrintError => 1}) 
	or die "Connection Error: $DBI::errstr\n";

#=============================================================================================
#=============================== Get options ===================================
#=============================================================================================
my $ProgramName  = "load_feat3.pl";
my ($temp,$project,$protein_id);

&GetOptions ("project:s" =>  \$project,"temp:s" =>  \$temp, "protein:s" => \$protein_id);

$temp = "seq_temp" if ($temp eq '');
my $tempname = $temp.".fsa";


unless (defined $project) { 
	print "#=========================================================\n";
	print "# ERROR: Project not defined, example: -project 1600genomes\n"; 
	print "#=========================================================\n";
	&usage;
}

unless (defined $protein_id) {
	print "#=========================================================\n";
	print "# ERROR: Protein ID not defined, example: -protein 555\n"; 
	print "#=========================================================\n";
}

sub usage {
	print "#=========================================================\n";
	print "# USAGE:\tperl $ProgramName -project <project name> -temp <name of temp FASTA file> -protein <protein id><\n";	
	print "# EXAMPLE:\tperl $ProgramName -project 1600genomes -temp temp_feat1\n";
	print "#=========================================================\n";
	exit( 1 );
}

#=============================================================================================
#=============================== Get list of protein IDs ===================================
#=============================================================================================

=put
my $prepare_protein_id = $dbh->prepare("SELECT protein_id from protein 
JOIN project_has_organism USING (org_id) 
WHERE projectname=?
ORDER BY RAND()");
WHERE projectname=? AND length < 1002 AND length > 1000");
$prepare_protein_id->execute($project)  or die $dbh->errstr;
$dbh->commit();
my @protein_ids;
while (my @row = $prepare_protein_id->fetchrow_array()) { push(@protein_ids, $row[0]) }
=cut

#=============================================================================================
#=============================== Print information to user ===================================
#=============================================================================================

print "#=========================================================\n";
print "# Project\t=\t $project\n";
print "# Protein ID\t=\t $protein_id\n";
print "#=========================================================\n";

#=============================================================================================
#=============================== Get protein sequence ===================================
#=============================================================================================
my ($seq, $length, $count) = ('','',0);

#----	SELECT sequence and length from protein 	----------------------------
my $prepare_seq = $dbh->prepare("SELECT seq,length from protein WHERE protein_id=?");
$prepare_seq->execute($protein_id)  or die $dbh->errstr;
$dbh->commit();

#=============================================================================================
#=============================== Analyze protein sequence ===================================
#=============================================================================================
if( my @row = $prepare_seq->fetchrow_array() ) { 
	$seq = $row[0];
	$seq =~ s/[ZXB\*]//g;
	$length = $row[1];

	#----	Check if protein has already been calculated	----------------------------
	#----	If protein WITH values is found (rows > 0), don't calculate again 
	my $prepare_duplicate = $dbh->prepare("SELECT * FROM feat WHERE protein_id=? AND seq_low_coverage IS NOT NULL AND signalp_p IS NOT NULL;");
	$prepare_duplicate->execute($protein_id) or die $dbh->err();
	$dbh->commit();

	#----	INSERT if not exists	----------------------------
	if ($prepare_duplicate->rows == 0) {

		#----	Temporary FASTA file	----------------------------
		open (TEMPFILE, ">", $tempname) or die "# ERROR: cannot create temp FASTA file, $tempname\n";
		print TEMPFILE ">protein_id:$protein_id\n$seq\n";

#=============================================================================================
#=============================== CALCULATE FEATURES ===================================
#=============================================================================================
		#----	Protparam	----------------------------
		my @protparam = `python $analysis_directory/protparam.py $tempname`;
		my ($protparam_values, $protparam_fields) = &parse_protparam(@protparam);
		#----	Psort	----------------------------
		my ($psort_n_fields, $psort_n_values,$psort_p_fields, $psort_p_values) = &psort();
		#----	Netphos	----------------------------
		my ($netphos_fields, $netphos_total,$netphos_serine, $netphos_threonine) = &netphos();
		#----	Signalp	----------------------------
		my ($signalp_n_fields,$signalp_n_values, $signalp_p_fields,$signalp_p_values) = &signalp();
		#----	tmhmm	----------------------------
		my ($tmhmm_fields, $tmhmm_values) = &tmhmm();
		#----	seg	----------------------------
		my ($seg_fields, $seg_low_total, $seg_low_avg, $seg_low_coverage, $seg_high_total, $seg_high_avg, $seg_high_coverage) = &seg($length);

#=============================================================================================
#=============================== MySQL INSERTS ===================================
#=============================================================================================
		my $all_feat = $protein_id.",".$protparam_values.",".$psort_n_values.",".$psort_p_values.",".$netphos_total.",".$netphos_serine.",".$netphos_threonine.",".$signalp_n_values.",".$signalp_p_values.",".$tmhmm_values.",".$seg_low_total.",".$seg_low_avg.",".$seg_low_coverage.",".$seg_high_total.",".$seg_high_avg.",".$seg_high_coverage;
		my $all_fields = "protein_id,".$protparam_fields.",".$psort_n_fields.",".$psort_p_fields.",".$netphos_fields.",".$signalp_n_fields.",".$signalp_p_fields.",".$tmhmm_fields.",".$seg_fields;

		#my $feat_count += () = $all_feat =~ /\,/g;
		#my $fields_count += () = $all_fields =~ /\,/g;
		#print "> $feat_count :: $fields_count\n";

		#print "$all_fields\n";
		#print "# $all_feat\n";

		my $query = "REPLACE INTO feat (".$all_fields.") VALUES (".$all_feat.");";
		my $prepare_feat = $dbh->prepare($query);
#		$prepare_feat->trace($dbh->parse_trace_flags('SQL|3|print'));
		$prepare_feat->execute() or die $dbh->err();
		$dbh->commit();
	} else {
		print "# SKIPPING: $protein_id\t";
	}
}
exit( 0 ); # Successful Execution

#=============================================================================================
#=============================== SUBRUTINE ===================================
#=============================================================================================

#--------------------------------------------------------------
# NetPhosBac
#--------------------------------------------------------------
sub netphos {
#--------------------------------------------------------------
	# Run program and error handling
	my @netphos = `$analysis_directory/ape $tempname | awk '\$6>0.5' `;
	my ($netphos_fields, $netphos_total,$netphos_serine, $netphos_threonine);
	if ($? == 0) {
		$netphos_total = scalar(@netphos);
		$netphos_serine = grep (/ S /, @netphos);
		$netphos_threonine = grep (/ T /, @netphos);
		$netphos_fields = 'netphos_total,netphos_serine,netphos_threonine';
	} else {
		($netphos_total,$netphos_serine,$netphos_threonine) = (0,0,0);
		$netphos_fields = 'netphos_total,netphos_serine,netphos_threonine';	
	}
	return ($netphos_fields, $netphos_total,$netphos_serine, $netphos_threonine);
}

#--------------------------------------------------------------
# Psort
#--------------------------------------------------------------
sub psort {
#--------------------------------------------------------------
	my($psort_n_fields, $psort_n_values, $psort_p_fields, $psort_p_values);
	#----------------
	# Psort	- Gram negatives
	#----------------
	# Run program and error handling
	my @psort_n = `perl $analysis_directory/ProtFun/bio-tools-psort-all/psort/bin/psort -n -o long $tempname `;	
	if ($? == 0) {
		# Basic processing
		my @psort_n_fields = split("\t",$psort_n[0]);		
		$psort_n[1] =~ s/\t$//g;		
		my @psort_n_values = split("\t",$psort_n[1]);		

		my (@psort_n_fields_short,@psort_n_values_short) = ((),());
		my @index_n = grep { $psort_n_fields[$_] =~ m/Score/ } 0..$#psort_n_fields;   		

		for my $i (@index_n) { 
			push (@psort_n_fields_short, $psort_n_fields[$i]."_n");
			push (@psort_n_values_short, $psort_n_values[$i]);
		}
		$psort_n_fields = join (",",@psort_n_fields_short);
		$psort_n_values = join (",",@psort_n_values_short);
	} else {
		&feature_error($protein_id, "PSORT N failed to run");
		$psort_n_fields = ("Cytoplasmic_Score_n,CytoplasmicMembrane_Score_n,Periplasmic_Score_n,OuterMembrane_Score_n,Extracellular_Score_n,Final_Score_n");
		$psort_n_values = ("0,0,0,0,0,0");
	}
	#----------------
	# Psort	- Gram positives
	#----------------
	# Run program and error handling
	my @psort_p = `perl $analysis_directory/ProtFun/bio-tools-psort-all/psort/bin/psort -p -o long $tempname `;
	if ($? == 0) {
		# Basic processing
		my @psort_p_fields = split("\t",$psort_p[0]);		
		$psort_p[1] =~ s/\t$//g;		
		my @psort_p_values = split("\t",$psort_p[1]);		

		my @index_p = grep { $psort_p_fields[$_] =~ m/Score/ } 0..$#psort_p_fields; 
	
		my (@psort_p_fields_short,@psort_p_values_short) = ((),());
		for my $i (@index_p) { 
			push (@psort_p_fields_short, $psort_p_fields[$i]."_p");
			push (@psort_p_values_short, $psort_p_values[$i]);
		}
		$psort_p_fields = join (",",@psort_p_fields_short);
		$psort_p_values = join (",",@psort_p_values_short);
	} else {
		&feature_error($protein_id, "PSORT P failed to run");
		$psort_p_fields = ("Cytoplasmic_Score_p,CytoplasmicMembrane_Score_p,Cellwall_Score_p,Extracellular_Score_p,Final_Score_p");
		$psort_p_values = ("0,0,0,0,0");
	}

	return($psort_n_fields, $psort_n_values, $psort_p_fields, $psort_p_values);
}

#--------------------------------------------------------------
# SignalP
#--------------------------------------------------------------
sub signalp {
#--------------------------------------------------------------
	my($signalp_n_fields,$signalp_n_values, $signalp_p_fields,$signalp_p_values);
	#----------------
	# SignalP - Gram positives	
	#----------------
	# Run program and error handling
	my @signalp_p = `$analysis_directory/signalp -t gram+ $tempname | grep -v "# SignalP"`;	

	if ($? == 0) {
		# Basic processing
		$signalp_p[0] =~ s/\#//g;		
		chomp $signalp_p[0];	
		chomp $signalp_p[1];

		my @signalp_p_fields = split(" ",$signalp_p[0]);	
		my @signalp_p_values = split(" ",$signalp_p[1]);	

		my @index_p = grep { $signalp_p_fields[$_] !~ m/(pos|name|Network|\?|Dmaxcut)/g } 0..$#signalp_p_fields;   		

		my (@signalp_p_fields_short,@signalp_p_values_short) = ((),());
		for my $i (@index_p) { 
			push (@signalp_p_fields_short, $signalp_p_fields[$i]."_p");
			push (@signalp_p_values_short, $signalp_p_values[$i]);
		}

		if ($signalp_p_values_short[4] >= 0.45) {
			push(@signalp_p_values_short,1);		
			push(@signalp_p_fields_short,"signalp_p");			
		} else {
			push(@signalp_p_values_short,0);		
			push(@signalp_p_fields_short,"signalp_p");			
		}
		$signalp_p_fields = join (",",@signalp_p_fields_short);
		$signalp_p_values = join (",",@signalp_p_values_short);
	} else {
		&feature_error($protein_id, "SignalP Gram Positive failed to run");
		$signalp_p_fields = ("Cmax_p,Ymax_p,Smax_p,Smean_p,D_p,signalp_p");
		$signalp_p_values = ("0,0,0,0,0,0");
	}

	#----------------
	# SignalP - Gram negatives	
	#----------------
	# Run program and error handling
	my @signalp_n = `$analysis_directory/signalp -t gram- $tempname | grep -v "# SignalP"`;
	
	if ($? == 0) {
		# Basic processing
		$signalp_n[0] =~ s/\#//g;		
		chomp $signalp_n[0];
		chomp $signalp_n[1];

		my @signalp_n_fields = split(" ",$signalp_n[0]);		
		my @signalp_n_values = split(" ",$signalp_n[1]);	

		my @index_n = grep { $signalp_n_fields[$_] !~ m/(pos|name|Network|\?|Dmaxcut)/g } 0..$#signalp_n_fields;   		

		my (@signalp_n_fields_short,@signalp_n_values_short) = ((),());

		for my $i (@index_n) { 
			push (@signalp_n_fields_short, $signalp_n_fields[$i]."_n");
			push (@signalp_n_values_short, $signalp_n_values[$i]);
		}
		if ($signalp_n_values_short[4] >= 0.45) {
			push(@signalp_n_values_short,1);		
			push(@signalp_n_fields_short,"signalp_n");			
		} else {
			push(@signalp_n_values_short,0);		
			push(@signalp_n_fields_short,"signalp_n");			
		}

		$signalp_n_fields = join (",",@signalp_n_fields_short);
		$signalp_n_values = join (",",@signalp_n_values_short);
	} else {
		&feature_error($protein_id, "SignalP Gram Negative failed to run");
		$signalp_n_fields = ("Cmax_n,Ymax_n,Smax_n,Smean_n,D_n,signalp_n");
		$signalp_n_values = ("0,0,0,0,0,0");
	}

	return($signalp_n_fields,$signalp_n_values, $signalp_p_fields,$signalp_p_values);
}

#--------------------------------------------------------------
# TMHMM
#--------------------------------------------------------------
sub tmhmm {
#--------------------------------------------------------------
	my ($tmhmm_fields, $tmhmm_values);
	# Run program and error handling
	my $tmhmm = `tmhmm -short $tempname`;	

	if ($? == 0) {
		# Basic processing
		chomp $tmhmm;
		$tmhmm =~ s/\:/=/g;

		my @tmhmm_sections = split("\t", $tmhmm);
		my (@tmhmm_fields,@tmhmm_values) = ((),());

		foreach my $section (@tmhmm_sections) {
			push(@tmhmm_fields,(split("=",$section))[0]);
			push(@tmhmm_values,(split("=",$section))[1]);
		}

		my @tmhmm_index = grep { $tmhmm_fields[$_] !~ m/(protein|len|Topology)/g } 0..$#tmhmm_fields;   		

		my (@tmhmm_fields_short,@tmhmm_values_short) = ((),());
		for my $i (@tmhmm_index) { 
			push (@tmhmm_fields_short, $tmhmm_fields[$i]);
			push (@tmhmm_values_short, $tmhmm_values[$i]);
		}
		$tmhmm_fields = join(",",@tmhmm_fields_short);
		$tmhmm_values = join(",",@tmhmm_values_short);	
	} else {
		&feature_error($protein_id, "TMHMM Failed to run");
		$tmhmm_fields = ("ExpAA,First60,PredHel");
		$tmhmm_values = ("0,0,0");	
	}
	return ($tmhmm_fields, $tmhmm_values);
}

#--------------------------------------------------------------
# SEG
#--------------------------------------------------------------
sub seg {
#--------------------------------------------------------------
	my $length = $_[0];
	my ($seg_fields, $seg_low_total, $seg_low_avg, $seg_low_coverage, $seg_high_total, $seg_high_avg, $seg_high_coverage);
	my @lengths = ();

	#----------------
	# SEG - high complexity
	#----------------
	# Run program and error handling
	my @seg_high = `$analysis_directory/seg $tempname -h | grep ">"`;	
	if ($? == 0) {
		$seg_high_total = scalar(@seg_high);
		foreach my $l (@seg_high) {
			my $t = (split(" ",$l))[0];
			$t =~ m/(\d+)-(\d+)/;
			push(@lengths,($2-$1));
		}
		$seg_high_avg = &average(\@lengths);
		$seg_high_coverage = ($seg_high_avg/$length)*100;
	} else {
		&feature_error($protein_id, "SEG High Failed to run");
		($seg_high_total, $seg_high_avg, $seg_high_coverage) = (0,0,0);
	}

	#----------------
	# SEG - low complexity
	#----------------
	# Run program and error handling
	@lengths = ();
	my @seg_low = `$analysis_directory/seg $tempname -l | grep ">"`;
	if ($? == 0) {
		$seg_low_total = scalar(@seg_low);
		foreach my $l (@seg_low) {
			my $t = (split(" ",$l))[0];
			$t =~ m/(\d+)-(\d+)/;
			push(@lengths,($2-$1));
		}
		$seg_low_avg = &average(\@lengths);
		$seg_low_coverage = ($seg_low_avg/$length)*100;
	} else {
		&feature_error($protein_id, "SEG Low has no hits or failed to run");
		($seg_low_total, $seg_low_avg, $seg_low_coverage) = (0,0,0);
	}
	$seg_fields ="seq_low_total, seq_low_avg, seq_low_coverage, seq_high_total, seq_high_avg, seq_high_coverage";

	return($seg_fields, $seg_low_total, $seg_low_avg, $seg_low_coverage, $seg_high_total, $seg_high_avg, $seg_high_coverage);
}

#--------------------------------------------------------------
# Protparam
#--------------------------------------------------------------
sub parse_protparam {
#--------------------------------------------------------------
	my (%features, $id, $v);
	my (@fields, @values, $fields_protparam, $values_protparam);
	foreach my $line (@_) {
		chomp $line;

		$line =~ s/:\s/:/g;	# Clean up the lines
		$line =~ s/\,//g;	# Clean up the lines
		next if $line =~ m/\#/;
		# Count the number of colons in line, more means multi value feature
		my $count = () = $line =~ m/:/gi; 
		
		# Get protein id
		if ($line =~ m/protein_id/) { $line =~ s/://g; $id = $line; next;	} 
		elsif ($line =~ m/^\n/) { $id = ''; next;	}

		if ($count > 1 and $line =~ /\{/) {	# If line is multi value feature
			my ($feat,$value) = split("{",$line);	# Get feature name and value
			
			$value =~ s/\{//g;	# Clean up lines
			$value =~ s/\}//g;	# Clean up lines
			$feat =~ s/ //g;	# Clean up lines
			$feat =~ s/  //g;	# Clean up lines
	
			my @array = split(" ", $value);	# Seperate sub values

			# Create new featurename = feature + subfeature and store each subvalue
			for my $a (@array) {
				my ($subfeat,$subvalue) = split(":", $a);
				$feat =~ s/:/_/g;
				$subfeat =~ s/\.//g;
				$features{$id}{$feat.$subfeat} = $subvalue;
			}		
		}	
		elsif ($line =~ m/:/) {	# If one-value feature
			my ($feat,$value) = split(":",$line);
			$feat =~ s/ //g;
			$feat =~ s/  //g;
			$features{$id}{$feat} = $value;
		}

		@fields = ( sort keys %{ $features{(sort keys %features)[0]} } );
		$fields_protparam = join(",", @fields);
	} 

	# Print values for each feature and protein
	foreach my $f ( sort keys %features ) {
		for my $r ( sort keys %{ $features{$f} } ) {
			if ($r =~ m/absorbtion/ or $r =~ m/extinction/ or $r =~ m/number/ or $r =~ m/total_atoms/) { 
				$features{$f}{$r} = ($features{$f}{$r})/($features{$f}{"sequence_length"});
			}
			$v =  sprintf "%0.3f", $features{$f}{$r};
			push (@values, $v);	
		}
	} 
	$values_protparam = join(",", @values);
	return ($values_protparam,$fields_protparam);
}
#--------------------------------------------------------------
# Array average
#--------------------------------------------------------------
sub average{
#--------------------------------------------------------------
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

#--------------------------------------------------------------
# Error reporting
#--------------------------------------------------------------
sub feature_error{
#--------------------------------------------------------------
	my ($pid, $message);
	($pid, $message) = @_;
	my $query = "INSERT INTO feature_errors (protein_id, message) VALUES ($pid, \"$message\")";
	my $prepare_feat = $dbh->prepare($query);
	$prepare_feat->execute();
	$dbh->commit();	
}
