#!/usr/bin/perl -w

#=============================================================================================
#=============================== USING STATEMENTS ============================================
#=============================================================================================

# STANDARD
use strict;
use warnings;

# BIOPERL SEQUENCE LOADING
use Bio::SeqIO;

# THREADING LIBRARY
use threads;
use threads::shared;
use Thread::Semaphore;

# COMMAND LINE ARGUMENTS
use Getopt::Long;

# DATA DUMPING
use Data::Dumper;

# MD5 HASHING
use Digest::MD5 qw(md5 md5_hex);

=put
#############################################################################
# Calculates protein features from protein sequence and stores as TAB file
#############################################################################

#--------------------------------------------------------------
# Input source
#--------------------------------------------------------------
Genefinder or protein sequecne download from database

#--------------------------------------------------------------
# Input
#--------------------------------------------------------------
head trial.fsa
>9832486_PF00004_23
VTSPQSAATDARNAMVAGMLASGISVNGLQPSHNPQVALQMFTTATTIDPSMCDAWLARVLAGDSDIEVLANAWATVRTFGWETRRLGLSDLEFHPEVSDGMFLRLRITSVESLAAAYAAALAEAKRYEEAAKLLDGIEPRNPFETELVGYVRGMLYFRTGRWPDVLKQFPAAAPWRQPELKAAAAAMATTALASLGVFEEATRRAQEAIEGDRVPGATNVALYTQGMCMRHLGREDEAAELLRRVYSRDAKFAPAREALDNLDYRLILTDPETIEARTDPWDPDSAPTREQAEAHRHAEEAARYLAEGDAELAAMLGMEQAKREIKLIKATTKVNLARTKMGLPVPVTSRHTLLLGPPGTGKTSVARAFTKQLCGLTVLRKPLVVETNRSKLLGRYMADAEKNTEELLEGALGGAVFFDEMHTLHERGYSQGDAYGNAIINTLLLYMENHRDELVVFGAGYANAMDKMLEVNQGLRRRFSTVIEFYSYTPPELLELTKMMGAENEDVVTDEAVEGLLPSYTKFYADENYSEDGDLIRGIDVLGNAGFVRNVVEKARDHRSFRLDDSELDAVLAGDVTEFSDEWLHKFKELTREDLAEGLSAAVAEKKPS
>9892486_PF00004_23
MTRPQAAAEDARNAMVAGLLASGISVNGLQPSHNPQVAAQMFTTATRLDPKMCDAWLARLLAGDQSIEVLAGAWAAVRTFGWETRRLGVTDLQFRPEVSDGLFLRLAITSVDSLACAYAAVLAEAKRYQEAAELLDATDPRHPFDAELVSYVRGVLYFRTKRWPDVLAQFPEATQWRHPELKAAGAAMATTALASLGVFEEAFRRAQEAIEGDRVPGAANIALYTQGMCLRHVGREEEAVELLRRVYSRDAKFTPAREALDNPNFRLILTDPETIEARTDPWDPDSAPTRAQTEAARHAEMAAKYLAEGDAELNAMLGMEQAKKEIKLIKSTTKVNLARAKMGLPVPVTSRHTLLLGPPGTGKTSVARAFTKQLCGLTVLRKPLVVETSRTKLLGRYMADAEKNTEEMLEGALGGAVFFDEMHTLHEKGYSQGDPYGNAIINTLLLYMENHRDELVVFGAGYAKAMEKMLEVNQGLRRRFSTVIEFFSYTPQELIALTQLMGRENEDVITEEESQVLLPSYTKFYMEQSYSEDGDLIRGIDLLGNAGFVRNVVEKARDHRSFRLDDEDLDAVLASDLTEFSEDQLRRFKELTREDLAEGLRAAVAEKKTK

#--------------------------------------------------------------
# Output
#--------------------------------------------------------------
head trial.fsa.tab
protein_id	protein_id	fasta	absorbtion_max	absorbtion_min	aliphatic_index	amino_acid_percentages_A	amino_acid_percentages_C	amino_acid_percentages_D	amino_acid_percentages_E	amino_acid_percentages_F	amino_acid_percentages_G	amino_acid_percentages_H	amino_acid_percentages_I	amino_acid_percentages_K	amino_acid_percentages_L	amino_acid_percentages_M	amino_acid_percentages_N	amino_acid_percentages_P	amino_acid_percentages_Q	amino_acid_percentages_R	amino_acid_percentages_S	amino_acid_percentages_T	amino_acid_percentages_V	amino_acid_percentages_W	amino_acid_percentages_Y	aromaticity	extinction_coeff_maxextinction_coeff_min	flexibility_average	gravy_coeff	in_vivo_halflife_ecoli	in_vivo_halflife_human	in_vivo_halflife_yeast	instability_index	isoelectric_point	molecular_weight	number_aromatic	number_negative	number_nonpolar	number_positive	number_uncharged	secondary_structure_helix	secondary_structure_sheet	secondary_structure_turn	sequence_length	total_atoms	total_molecules	weight_aromatic	weight_negative	weight_nonpolar	weight_positive	weight_uncharged	Cytoplasmic_Score_n	CytoplasmicMembrane_Score_n	Periplasmic_Score_n	OuterMembrane_Score_n	Extracellular_Score_n	Final_Score_n	Cytoplasmic_Score_p	CytoplasmicMembrane_Score_p	Cellwall_Score_p	Extracellular_Score_p	Final_Score_p	netphos_total	netphos_serine	netphos_threonine	Cmax_n	Ymax_n	Smax_n	Smean_n	D_n	signalp_n	Cmax_p	Ymax_p	Smax_p	Smean_p	D_p	signalp_p	ExpAA	First60	PredHel	seq_low_total	seq_low_avg	seq_low_coverage	seq_high_total	seq_high_avg	seq_high_coverage
9832486_PF00004_23	9832486_PF00004_23	trial.fsa	0.002	0.002	84.574	0.134	0.005	0.061	0.093	0.033	0.067	0.018	0.025	0.038	0.111	0.031	0.031	0.043	0.021	0.074	0.044	0.067	0.062	0.011	0.030	0.074	107.389	107.082	1.001	-0.285	600	6000	1200	45.160	5.008	67286.210	0.074	0.154	0.523	0.130	0.193	0.272	0.370	0.185	610	18.428	610	7994.830	13311.110	36648.150	12908.130	15393	8.96	0.51	0.26	0.01	0.26	8.96	7.501.15	0.62	0.73	7.50	1079	390	625	0.115	0.175	0.414	0.275	0.222	0	0.142	0.221	0.501	0.290	0.248	0	0.60	0.59	0	3	16.6666666666667	2.73224043715847	4	138.25	22.6639344262295


#--------------------------------------------------------------
# USAGE
#--------------------------------------------------------------
perl CMGfunc_calcFeatures.pl -fasta genomeA.proteins.fsa

Creates file.fsa.tab
=cut

#=============================================================================================
#=============================== Set defaults ================================================
#=============================================================================================
my $ProtFun_path = `locate -b ProtFun | grep -v -P "ProtFun."`; chomp $ProtFun_path;

#=============================================================================================
#=============================== Get options =================================================
#=============================================================================================

my $fasta;
my $i_counter = 0;

&GetOptions ("fasta:s" =>  \$fasta);

unless (defined $fasta) { 
	print "#==================================================================================\n";
	printf ("%-40s\t==>\t%-50s\n", "# ERROR", "Protein FASTA file not defined, example: -fasta genomeA.proteins.fsa"); 
	&usage;
}

sub usage {
	print "#==================================================================================\n";
	printf ("%-40s\t==>\t%-50s\n", "# USAGE", "perl CMGfunc_calcFeatures.pl -fasta <name of FASTA file>");	
	printf ("%-40s\t==>\t%-50s\n", "# EXAMPLE", "perl CMGfunc_calcFeatures.pl -fasta genomeA.proteins.fsa");
	print "#==================================================================================\n";
	exit( 1 );
}

#=============================================================================================
#=============================== Definitions for Multithread =================================
#=============================================================================================

# CPU Information -- Sorry for this hack, but I wanted to avoid adding another perl module
# lscpu shows cpu information
#  -> grep for line showing number of cpus
#      -> reverse the result
#          -> take the first field
#              -> reverse the result back to normal
# (Or 1 if the command for some reason comes up empty)
my $n_threads = scalar(`lscpu | grep 'CPU(s):' | rev | cut -d " " -f 1 | rev`) || 1;

# Queue Sizes
my $calc_queue_max  = $n_threads * 5;          # Maximum number of slots in the feature processing queue
my $print_queue_max = $n_threads * 5;          # Maximum number of slots in the printing queue

# Calculation Queue
my $calc_queue_free = Thread::Semaphore->new($calc_queue_max);
my $calc_queue_used = Thread::Semaphore->new(0);
my @calc_queue :shared; # Important Note: Queued info is three consecutive items!

# Printing Queue
my $print_queue_free = Thread::Semaphore->new($print_queue_max);
my $print_queue_used = Thread::Semaphore->new(0);
my @print_queue :shared;

# Thread States
my $load_thread_state :shared = 0; # (0=Running,1=Done)
my $calc_thread_state = Thread::Semaphore->new((-1 * $n_threads) + 1);

#=============================================================================================
#=============================== Subroutines for Multithread =================================
#=============================================================================================
sub loadSequences{
	# Get file from input argument
	my $fname = shift || die("sub loadSequences requires an argument to run on");
	# Open Bio::SeqIO object for parsing file
	my $seqIO_obj = Bio::SeqIO->new(-file=>"<$fname");
	# For each sequence in the file
	while( my $seq_obj = $seqIO_obj->next_seq ) {
		# Wait until there is space in the calculation queue
		$calc_queue_free->down();
		{
			lock(@calc_queue);
			# Write the sequence to the calculation queue - Again, note three items added
			push(@calc_queue, ($seq_obj->id, $seq_obj->length, $seq_obj->seq()) );
		}
		# Inform feature processing threads that there is new info to process
		$calc_queue_used->up();
	}
	# Inform feature processing threads that there will be no more info from you
	lock($load_thread_state);
	$load_thread_state = 1;
}

sub calcSequences{
	# Done loading features?
	my $load_done = 0;
	# Get current thread ID and name the tempfile accordingly
	my $thread_id = threads->tid();
	my $tempname = "temp_$thread_id.fasta";
	# Create tempfile (as a check, and so we clean it later)
	open (TEMPFILE, ">", $tempname) or die printf ("%-40s\t==>\t%-50s\n", "# ERROR", "cannot create temp FASTA file, $tempname");
	# While true...
	while( 1 ){
		# Get whether the loading thread is done submitting info
		{
			lock($load_thread_state);
			$load_done = $load_thread_state;
		}
		# If there is info in the processing queue
		if( $calc_queue_used->down_nb() ) {
			# Get the info for processing
			my ($protein_id, $length, $seq);
			{
				lock(@calc_queue);
				$protein_id = shift(@calc_queue);
				$length = shift(@calc_queue);
				$seq = shift(@calc_queue);
			}
			# Tell the loader there is space again
			$calc_queue_free->up();

			# Write the temporary file for this thread
			seek(TEMPFILE,0,0);
			truncate(TEMPFILE,0);
			$protein_id =~ s/\_/ /g; # Some programs crash if the name is to long!
			$seq =~ s/[ZXB\*]//g;
			$seq =~ s/U/C/g;
			print TEMPFILE ">protein_id:$protein_id\n$seq\n";

        		$i_counter++;
			print "# Running sequence number: $i_counter\n" if ($i_counter % 10 == 0);

#=============================================================================================
#=============================== CALCULATE FEATURES ==========================================
#=============================================================================================
			#.............. Protparam ..............
			my @protparam = `python $ProtFun_path/protparam.py $tempname`;
			my ($protparam_values, $protparam_fields) = &parse_protparam(@protparam);

			#.............. Psort ..............
			my ($psort_n_fields, $psort_n_values,$psort_p_fields, $psort_p_values) = &psort($tempname);

			#.............. Netphos ..............
			my ($netphos_fields, $netphos_total,$netphos_serine, $netphos_threonine) = &netphos($tempname);

			#.............. Signalp ..............
			my ($signalp_n_fields,$signalp_n_values, $signalp_p_fields,$signalp_p_values) = &signalp($tempname);

			#.............. tmhmm ..............
			my ($tmhmm_fields, $tmhmm_values) = &tmhmm($tempname);

			#.............. seg ..............
			my ($seg_fields, $seg_low_total, $seg_low_avg, $seg_low_coverage, $seg_high_total, $seg_high_avg, $seg_high_coverage) = &seg($tempname, $length);
#=============================================================================================
#=============================== FEATURE VECTOR ==============================================
#=============================================================================================
			# Build concatanated string with commas
			my $all_feat = $protein_id.",".$protein_id.",".$fasta.",".$protparam_values.",".$psort_n_values.",".$psort_p_values.",".$netphos_total.",".$netphos_serine.",".$netphos_threonine.",".$signalp_n_values.",".$signalp_p_values.",".$tmhmm_values.",".$seg_low_total.",".$seg_low_avg.",".$seg_low_coverage.",".$seg_high_total.",".$seg_high_avg.",".$seg_high_coverage;
			my $all_fields = "protein_id,"."protein_nr,"."fasta,".$protparam_fields.",".$psort_n_fields.",".$psort_p_fields.",".$netphos_fields.",".$signalp_n_fields.",".$signalp_p_fields.",".$tmhmm_fields.",".$seg_fields;

			# Format data
			$all_feat =~ s/\,/\t/g;
			$all_feat =~ s/\.000//g;
			$all_fields =~ s/\,/\t/g;

			# Wait for space in the printing queue
			$print_queue_free->down();
			# Add the features to the writing queue
			{
				lock(@print_queue);
				push(@print_queue, $all_fields);
				push(@print_queue, $all_feat);
			}
			# Tell the printing thread there is new info in the queue
			$print_queue_used->up();
		} elsif( $load_done ) {
			# There was nothing in the queue and the loading thread was done
			# Tell the writer you're done, clean up and terminate
			$calc_thread_state->up();
			unlink $tempname;
			return;
		} else {
			# Nothing in the queue but still loading info
			# Yield processing time to other threads
			threads->yield();
		}
	}
}

sub printSequences{
	# Output file name from argument
	my $fname = shift || die("sub loadSequences requires an argument to run on");
	# Open output file
	open (OUTFILE, ">", $fname.".tab");
	# Whether the calculation threads are done
	my $calc_done = 0;
	# First Print Call?
	my $is_first_print = 1;
	# While true...
	while( 1 ){
		# Check whether calculation threads are done
		if( $calc_thread_state->down_nb() ) {
			$calc_done = 1;
		}
		# Check whether there is data to print
		if( $print_queue_used->down_nb() ) {
			# Get the data to print
			my $print_fields;
			my $print_features;
			{
				lock(@print_queue);
				$print_fields = shift(@print_queue);
				$print_features = shift(@print_queue);
			}
			# Tell the calculation threads there is space in the printing queue
			$print_queue_free->up();
			if( $is_first_print ) {
				print OUTFILE $print_fields, "\n";
				$is_first_print = 0;
			}
			print OUTFILE $print_features, "\n";
		} elsif( $calc_done ) {
			# No data to print and calculation threads are all done
			# Therefore, terminate
			return;
		} else {
			# No data to print, but calculation threads are running
			# Yield processing time to other threads
			threads->yield();
		}
	}	
}

#=============================================================================================
#=============================== THREAD CREATION =============================================
#=============================================================================================

my $loading_thread = threads->create(\&loadSequences, $fasta);
my @calc_threads = (0 .. ($n_threads-1));
for( my $i = 0; $i < $n_threads; $i++ ) {
	$calc_threads[$i] = threads->create(\&calcSequences);
}
printSequences($fasta);

# Make sure threads all rejoin
$loading_thread->join();
foreach my $t (@calc_threads)
{
	$t->join();
}

#=============================================================================================
#=============================== FEATURE SUBROUTINES =========================================
#=============================================================================================

#......................................................................
#.............. NetPhosBac ..............
#......................................................................
sub netphos {
	my $tempname = shift;
	#.............. Run program and error handling ..............
	my @netphos = `$ProtFun_path/ape $tempname | awk '\$6>0.5' `;
	my ($netphos_fields, $netphos_total,$netphos_serine, $netphos_threonine);
	if ($? == 0) {
		$netphos_total = scalar(@netphos);
		$netphos_serine = grep (/ S /, @netphos);
		$netphos_threonine = grep (/ T /, @netphos);
		$netphos_fields = 'netphos_total,netphos_serine,netphos_threonine';
	} 
	else {
		($netphos_total,$netphos_serine,$netphos_threonine) = ("0,0,0");
		$netphos_fields = 'netphos_total,netphos_serine,netphos_threonine';	
	}
	return ($netphos_fields, $netphos_total,$netphos_serine, $netphos_threonine);
}

#......................................................................
#.............. Psort..............
#......................................................................
sub psort {
	my $tempname = shift;
	my($psort_n_fields, $psort_n_values, $psort_p_fields, $psort_p_values);
	#.............. Psort - Gram negatives ..............
	#.............. Run program and error handling ..............

	my @psort_n = `perl $ProtFun_path/ProtFun/bio-tools-psort-all/psort/bin/psort -n -o long $tempname > /dev/null 2>&1`;	
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
	} 
	else {
		printf ("%-40s\t==>\t%-50s\n", "# WARNING", "PSORT N failed to run");
		$psort_n_fields = ("Cytoplasmic_Score_n,CytoplasmicMembrane_Score_n,Periplasmic_Score_n,OuterMembrane_Score_n,Extracellular_Score_n,Final_Score_n");
		$psort_n_values = ("0,0,0,0,0,0");
	}
	#..............Psort - Gram positives ..............
	#.............. Run program and error handling ..............

	my @psort_p = `perl $ProtFun_path/ProtFun/bio-tools-psort-all/psort/bin/psort -p -o long $tempname > /dev/null 2>&1`;
	if ($? == 0) {
		#.............. Basic processing ..............
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
	} 
	else {
		printf ("%-40s\t==>\t%-50s\n", "# WARNING", "PSORT P failed to run");
		$psort_p_fields = ("Cytoplasmic_Score_p,CytoplasmicMembrane_Score_p,Cellwall_Score_p,Extracellular_Score_p,Final_Score_p");
		$psort_p_values = ("0,0,0,0,0");
	}

	return($psort_n_fields, $psort_n_values, $psort_p_fields, $psort_p_values);
}

#......................................................................
#.............. SignalP ..............
#......................................................................
sub signalp {
	my $tempname = shift;
	my($signalp_n_fields,$signalp_n_values, $signalp_p_fields,$signalp_p_values);
	#.............. SignalP - Gram positives..............
	#.............. Run program and error handling ..............

	my @signalp_p = `$ProtFun_path/signalp -t gram+ $tempname | grep -v "# SignalP"`;	

	if ($? == 0) {
		#.............. Basic processing ..............
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
		} 
		else {
			push(@signalp_p_values_short,0);		
			push(@signalp_p_fields_short,"signalp_p");			
		}
		$signalp_p_fields = join (",",@signalp_p_fields_short);
		$signalp_p_values = join (",",@signalp_p_values_short);
	} 
	else {
		printf ("%-40s\t==>\t%-50s\n", "# WARNING", "SignalP Gram Positive failed to run");
		$signalp_p_fields = ("Cmax_p,Ymax_p,Smax_p,Smean_p,D_p,signalp_p");
		$signalp_p_values = ("0,0,0,0,0,0");
	}

	#.............. SignalP - Gram negatives..............
	#.............. Run program and error handling ..............
	my @signalp_n = `$ProtFun_path/signalp -t gram- $tempname | grep -v "# SignalP"`;
	
	if ($? == 0) {
		#.............. Basic processing ..............
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
		} 
		else {
			push(@signalp_n_values_short,0);		
			push(@signalp_n_fields_short,"signalp_n");			
		}

		$signalp_n_fields = join (",",@signalp_n_fields_short);
		$signalp_n_values = join (",",@signalp_n_values_short);
	} 
	else {
		printf ("%-40s\t==>\t%-50s\n", "# WARNING", "SignalP Gram Negative failed to run");
		$signalp_n_fields = ("Cmax_n,Ymax_n,Smax_n,Smean_n,D_n,signalp_n");
		$signalp_n_values = ("0,0,0,0,0,0");
	}

	return($signalp_n_fields,$signalp_n_values, $signalp_p_fields,$signalp_p_values);
}

#......................................................................
#.............. TMHMM ..............
#......................................................................
sub tmhmm {
	my $tempname = shift;
	my ($tmhmm_fields, $tmhmm_values);
	#.............. Run program and error handling ..............
	my $tmhmm = `$ProtFun_path/tmhmm-2.0c/bin/tmhmm -short $tempname`;	

	if ($? == 0 or $tmhmm == '') {	# $? == 0, this error handling is not working. When tmhmm returns nothing, the else loop is not executed.
		#.............. Basic processing ..............
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
		printf ("%-40s\t==>\t%-50s\n", "# WARNING", "TMHMM Failed to run");
		$tmhmm_fields = ("ExpAA,First60,PredHel");
		$tmhmm_values = ("0,0,0");	
	}
	
	return ($tmhmm_fields, $tmhmm_values);
}

#......................................................................
#.............. SEG ..............
#......................................................................
sub seg {
	# NOTE: params - tempname, length
	my $tempname = shift;
	my $length   = shift;
	my ($seg_fields, $seg_low_total, $seg_low_avg, $seg_low_coverage, $seg_high_total, $seg_high_avg, $seg_high_coverage);
	my @lengths = ();

	#..............SEG - high complexity ..............
	#.............. Run program and error handling ..............
	my @seg_high = `$ProtFun_path/seg $tempname -h | grep ">"`;	
	if ($? == 0) {
		$seg_high_total = scalar(@seg_high);
		foreach my $l (@seg_high) {
			my $t = (split(" ",$l))[0];
			$t =~ m/\((\d+)-(\d+)\)/;
			push(@lengths,($2-$1));
		}
		$seg_high_avg = &average(\@lengths);
		$seg_high_coverage = ($seg_high_avg/$length)*100;

	} else {
		printf ("%-40s\t==>\t%-50s\n", "# WARNING", "SEG High Failed to run");
		($seg_high_total, $seg_high_avg, $seg_high_coverage) = ("0","0","0");
	}

	#.............. SEG - low complexity ..............
	#.............. Run program and error handling ..............
	@lengths = ();
	my @seg_low = `$ProtFun_path/seg $tempname -l | grep ">"`;
	if ($? == 0) {
		$seg_low_total = scalar(@seg_low);
		foreach my $l (@seg_low) {
			my $t = (split(" ",$l))[0];
			$t =~ m/\((\d+)-(\d+)\)/;
			push(@lengths,($2-$1));
		}
		$seg_low_avg = &average(\@lengths);
		$seg_low_coverage = ($seg_low_avg/$length)*100;
	} else {
		#printf ("%-40s\t==>\t%-50s\n", "# WARNING", "SEG Low has no hits or failed to run");
		($seg_low_total, $seg_low_avg, $seg_low_coverage) = ("0","0","0");
	}
	$seg_fields ="seq_low_total,seq_low_avg,seq_low_coverage,seq_high_total,seq_high_avg,seq_high_coverage";

	return($seg_fields, $seg_low_total, $seg_low_avg, $seg_low_coverage, $seg_high_total, $seg_high_avg, $seg_high_coverage);
}

#......................................................................
#.............. Protparam ..............
#......................................................................
sub parse_protparam {

	my (%features, $id, $v);
	my (@fields, @values, $fields_protparam, $values_protparam);
	foreach my $line (@_) {
		chomp $line;

		$line =~ s/:\s/:/g;	# Clean up the lines
		$line =~ s/\,//g;	# Clean up the lines
		next if $line =~ m/\#/;
		#.............. Count the number of colons in line, more means multi value feature ..............
		my $count = () = $line =~ m/:/gi; 
		
		#.............. Get protein id ..............
		if ($line =~ m/protein_id/) { $line =~ s/://g; $id = $line; next;	} 
		elsif ($line =~ m/^\n/) { $id = ''; next;	}

		if ($count > 1 and $line =~ /\{/) {	# If line is multi value feature
			my ($feat,$value) = split("{",$line);	# Get feature name and value
			
			$value =~ s/\{//g;	# Clean up lines
			$value =~ s/\}//g;	# Clean up lines
			$feat =~ s/ //g;	# Clean up lines
			$feat =~ s/  //g;	# Clean up lines
	
			my @array = split(" ", $value);	# Seperate sub values

			#.............. Create new featurename = feature + subfeature and store each subvalue ..............
			for my $a (@array) {
				my ($subfeat,$subvalue) = split(":", $a);
				$feat =~ s/:/_/g;
				$subfeat =~ s/\.//g;
				$features{$id}{$feat.$subfeat} = $subvalue;
			}		
		}	
		elsif ($line =~ m/:/) {	# If one-value feature
			my ($feat,$value) = split(":",$line);
			$feat =~ s/ //g; $feat =~ s/  //g;
			$features{$id}{$feat} = $value;
		}

		@fields = ( sort keys %{ $features{(sort keys %features)[0]} } );
		$fields_protparam = join(",", @fields);
	} 

	#.............. Print values for each feature and protein ..............
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
#......................................................................
#.............. Array average ..............
#......................................................................
sub average{

        my($data) = @_;
        if (not @$data)  {	die("Empty array\n")	}
        my $total = 0;
        foreach (@$data) {	$total += $_		}
        my $average = $total / @$data;
        return $average;
}
