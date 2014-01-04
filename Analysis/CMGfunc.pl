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
Genefinder or protein sequecne download from database

#--------------------------------------------------------------
# Input
#--------------------------------------------------------------
Protein FASTA file
head TRIALDATA/trial.fsa
>9832486_PF00004_23
VTSPQSAATDARNAMVAGMLASGISVNGLQPSHNPQVALQMFTTATTIDPSMCDAWLARVLAGDSDIEVLANAWATVRTFGWETRRLGLSDLEFHPEVSDGMFLRLRITSVESLAAAYAAALAEAKRYEEAAKLLDGIEPRNPFETELVGYVRGMLYFRTGRWPDVLKQFPAAAPWRQPELKAAAAAMATTALASLGVFEEATRRAQEAIEGDRVPGATNVALYTQGMCMRHLGREDEAAELLRRVYSRDAKFAPAREALDNLDYRLILTDPETIEARTDPWDPDSAPTREQAEAHRHAEEAARYLAEGDAELAAMLGMEQAKREIKLIKATTKVNLARTKMGLPVPVTSRHTLLLGPPGTGKTSVARAFTKQLCGLTVLRKPLVVETNRSKLLGRYMADAEKNTEELLEGALGGAVFFDEMHTLHERGYSQGDAYGNAIINTLLLYMENHRDELVVFGAGYANAMDKMLEVNQGLRRRFSTVIEFYSYTPPELLELTKMMGAENEDVVTDEAVEGLLPSYTKFYADENYSEDGDLIRGIDVLGNAGFVRNVVEKARDHRSFRLDDSELDAVLAGDVTEFSDEWLHKFKELTREDLAEGLSAAVAEKKPS
>9892486_PF00004_23
MTRPQAAAEDARNAMVAGLLASGISVNGLQPSHNPQVAAQMFTTATRLDPKMCDAWLARLLAGDQSIEVLAGAWAAVRTFGWETRRLGVTDLQFRPEVSDGLFLRLAITSVDSLACAYAAVLAEAKRYQEAAELLDATDPRHPFDAELVSYVRGVLYFRTKRWPDVLAQFPEATQWRHPELKAAGAAMATTALASLGVFEEAFRRAQEAIEGDRVPGAANIALYTQGMCLRHVGREEEAVELLRRVYSRDAKFTPAREALDNPNFRLILTDPETIEARTDPWDPDSAPTRAQTEAARHAEMAAKYLAEGDAELNAMLGMEQAKKEIKLIKSTTKVNLARAKMGLPVPVTSRHTLLLGPPGTGKTSVARAFTKQLCGLTVLRKPLVVETSRTKLLGRYMADAEKNTEEMLEGALGGAVFFDEMHTLHEKGYSQGDPYGNAIINTLLLYMENHRDELVVFGAGYAKAMEKMLEVNQGLRRRFSTVIEFFSYTPQELIALTQLMGRENEDVITEEESQVLLPSYTKFYMEQSYSEDGDLIRGIDLLGNAGFVRNVVEKARDHRSFRLDDEDLDAVLASDLTEFSEDQLRRFKELTREDLAEGLRAAVAEKKTK

#--------------------------------------------------------------
# Output
#--------------------------------------------------------------
head TRIALDATA/trial.fsa.tab.vec.res 
Testseq:	9892486_PF00004_23 	Network:	proteindata_clan_135.tab.A.net 	Score:	0.830191098163
Testseq:	9892486_PF00004_23 	Network:	proteindata_clan_174.tab.A.net 	Score:	0.805813615358

#--------------------------------------------------------------
# USAGE
#--------------------------------------------------------------
perl CMGfunc.pl -fasta file.fsa -net /home/cmgfunc/CMGfunc/NETS -mse 0.8

=cut

#=======================================================================================
#=======================================================================================
#==============================	ERROR HANDLING	========================================
#=======================================================================================
#=======================================================================================

#.............. Inputs ..............
my $CMGfunc_path	= `locate CMGfunc | grep -P "CMGfunc/indirect" | grep -P -v "CMGfunc/indirect." | grep -v MySQL`; chomp $CMGfunc_path;
my $pybrain_path	= `locate -b pybrain | grep -v "pybrain." | grep -v MySQL | grep -v local`; chomp $pybrain_path;
my $ProtFun_path	= `locate -b ProtFun | grep -v -P "ProtFun."`; chomp $ProtFun_path;
my $defaultNet_path	= `locate -b NETS | grep -v MySQL`; chomp $defaultNet_path;
unless (-d $defaultNet_path) { printf ("%-40s:\t%-50s\n", "# ERROR checkPATHS", "defaultNet_directory : $defaultNet_path");	exit }

#.............. Get options ..............
my ($fasta, $NET_path, $set_score, $clean);
&GetOptions ("fasta:s" =>  \$fasta, "net:s" => \$NET_path, "clean:f" => \$clean, "score:f" => \$set_score);

#.............. Input files and paths exists ..............
if (!defined $NET_path)	{	
	printf ("%-30s:\t%-50s\n", "# WARNING", "Path to networkfiles not defined, default is used: $defaultNet_path");			
	$NET_path = $defaultNet_path."/"; 
}
if (!defined $fasta)	{ printf ("%-30s:\t%-50s\n", "# ERROR", "Protein FASTA file not defined, example: perl CMGfunc.pl -fasta genomeA.proteins.fsa");	&usage;	exit; }
if (!defined $set_score){ $set_score	= 0.8;	}
if (!defined $clean)	{ $clean		= 1;	}

#.............. Test that paths exists or exit ..............
&checkPATHS();

#.............. Test that scripts exists in the right path or exit ..............
&checkSCRIPTS();

#.............. Print usage ..............
sub usage {	
	printf ("%-30s:\t%-50s\n", "# USAGE", "perl CMGfunc.pl -fasta <name of FASTA file> -net <path to network files> -mse <MSE cutoff for raw network results>");	
	printf ("%-30s:\t%-50s\n", "# EXAMPLE", "perl CMGfunc.pl -fasta genomeA.proteins.fsa -net /home/cmgfunc/CMGfunc/NETS -mse 0.01");
	printf ("%-30s:\t%-50s\n", "# INFO", "mse is optional, default is 0.04, cutoff is used to filter raw network results");
	exit( 1 );	}


#.............. Test that file is protein FASTA or exit ..............
unless (&checkFileTypeFASTA($fasta) eq "FASTA")		{	printf ("%-30s:\t%-50s\n", "# ERROR", "File, $fasta, is not FASTA format");	exit }

#=======================================================================================
#=======================================================================================
#==================================	INFORMATION	========================================
#=======================================================================================
#=======================================================================================
my $fsa_seq_count	= `grep -c '>' $fasta`; chomp $fsa_seq_count;
 
#.............. Print information to user ..............
print "#==================================================================================\n";
printf ("%-30s:\t%-60s|\n", '# Protein FASTA file', 		$fasta);
printf ("%-30s:\t%-60s|\n", '# Seqs. in FASTA file', 		$fsa_seq_count); 
printf ("%-30s:\t%-60s|\n", '# Network Score cutoff', 		$set_score);
printf ("%-30s:\t%-60s|\n", '# Python pybrain library path',$pybrain_path);
printf ("%-30s:\t%-60s|\n", '# Analysis directory path', 	$ProtFun_path);
printf ("%-30s:\t%-60s|\n", '# Path to network files', 		$NET_path);
print "#==================================================================================\n";

#=======================================================================================
#=======================================================================================
#=======================================	MAIN	================================
#=======================================================================================
#=======================================================================================

if ($clean == 0) {
	printf ("%-30s\t==>\t%-50s\n", "# OPTION -clean defined", "REMOVING files $fasta.tab $fasta.tab.vec $fasta.tab.vec.res $fasta.tab.vec.res.func");
	`rm -f $fasta.tab $fasta.tab.vec $fasta.tab.vec.res $fasta.tab.vec.res.func`;
}

#.............. TEST Calculate features output ..............
if (-e $fasta.".tab" and &checkFileTypeTAB($fasta, $fasta.".tab") eq "TRUE") {
	printf ("%-40s\t==>\t%-50s\n", "# FILE", "already exists and is the right format: $fasta.tab");

	&normFeatures($fasta);		#.............. Normalize features ..............
	unless (&checkFileTypeVECTOR($fasta, $fasta.".tab.vec") eq "TRUE") 					{ exit }

	&testNetworks($fasta);		#.............. Compare networks ..............
	unless (&checkFileTypeRESULTS($fasta, $fasta.".tab.vec.res", $NET_path) eq "TRUE")	{ exit }
} 

#.............. In the case that the first file is not right, do everything ..............
else {
	
	&calcFeatures($fasta); 		#.............. Calculate features ..............
	unless (&checkFileTypeTAB($fasta, $fasta.".tab") eq "TRUE") 						{ exit }

	&normFeatures($fasta);		#.............. Normalize features ..............
	unless (&checkFileTypeVECTOR($fasta, $fasta.".tab.vec") eq "TRUE")					{ exit }

	&testNetworks($fasta);		#.............. Compare networks ..............
	unless (&checkFileTypeRESULTS($fasta, $fasta.".tab.vec.res", $NET_path) eq "TRUE")	{ exit }
}

#=======================================================================================
#=======================================================================================
#==========================	CALCULATION SUBRUTINES	====================================
#=======================================================================================
#=======================================================================================

#.............. calcFeatures ..............
sub calcFeatures {
	my $fasta = $_[0];
	printf ("%-40s\t==>\t%-50s\n", "# RUNNING CMGfunc_calcFeatures.pl", "will create file: $fasta.tab");
	system ("perl $CMGfunc_path/CMGfunc_calcFeatures.pl -fasta $fasta");
	unless ($? == 0)			{ printf ("%-40s\t==>\t%-50s\n", "# ERROR", "$CMGfunc_path/CMGfunc_calcFeatures.pl did not run with success");		exit }
	unless (-e $fasta.".tab")	{ printf ("%-40s\t==>\t%-50s\n", "# ERROR", "$CMGfunc_path/CMGfunc_calcFeatures.pl did not return file $fasta.tab");	exit }
}

#.............. normFeatures ..............
sub normFeatures {
	printf ("%-40s\t==>\t%-50s\n", "# RUNNING CMGfunc_normFeatures.pl", "will create file: $fasta.tab.vec");
	system ("perl $CMGfunc_path/CMGfunc_normFeatures.pl -tab $fasta.tab");
	unless ($? == 0)				{ printf ("%-40s\t==>\t%-50s\n", "# ERROR", "$CMGfunc_path/CMGfunc_normFeatures.pl did not run with success");				exit }
	unless (-e $fasta.".tab.vec")	{ printf ("%-40s\t==>\t%-50s\n", "# ERROR", "$CMGfunc_path/CMGfunc_normFeatures.pl did not return file $fasta.tab.vec");	exit }
}

#.............. testNetworks ..............
sub testNetworks {
	printf ("%-40s\t==>\t%-50s\n", "# RUNNING CMGfunc_testNetworks.py", "will create file: $fasta.tab.vec.res");
	printf ("%-40s\t==>\t%-50s\n", "# RUNNING CMGfunc_testNetworks.py", "uses $NET_path");
	system ("python $CMGfunc_path/CMGfunc_testNetworks.py -i $fasta.tab.vec -n $NET_path -m $set_score ");
	system ("sort -k1 $fasta.tab.vec.res > tmp");
	system ("mv tmp $fasta.tab.vec.res");
	unless ($? == 0)					{ printf ("%-40s\t==>\t%-50s\n", "# ERROR", "$CMGfunc_path/CMGfunc_testNetworks.py did not run with success");					exit }
	unless (-e $fasta.".tab.vec.res")	{ printf ("%-40s\t==>\t%-50s\n", "# ERROR", "$CMGfunc_path/CMGfunc_testNetworks.py did not return file $fasta.tab.vec.res");	exit }
}
#=======================================================================================
#=======================================================================================
#===========================	FILETYPE SUBRUTINES	====================================
#=======================================================================================
#=======================================================================================

#.............. checkFileTypeRESULTS ..............
sub checkFileTypeRESULTS {
	my $status			= "TRUE";
	my $res_filename	= $_[1];
	my $fsa_filename	= $_[0];
	my $net_library		= $_[2];
	my $i				= 0;

	# Open the TAB file and test each vector for problems
	open(FP0,$res_filename);
	while( my $line = <FP0> )  {
		next if $line =~ m/^#/ ;
		my @elements = split("\t", $line);

		unless (scalar(@elements) == 6)	{ 
			printf ("%-40s\t==>\t%-50s\n", "# ERROR checkFileTypeRESULTS", "RES vectore, $res_filename, in line $i does not have the right number of elements");
			$status = "FALSE"; return $status; }

		foreach my $item (@elements)		{ 
			if ($item eq '') { printf ("%-40s\t==>\t%-50s\n", "# ERROR checkFileTypeRESULTS", "RES vectore, $res_filename, in line $i has empty elements");		
			$status = "FALSE"; return $status; }}
		$i++;
	}
	close(FP0);
	return $status;
}

#.............. checkFileTypeVECTOR ..............
sub checkFileTypeVECTOR {
	my $status				= "TRUE";
	my $vector_filename 	= $_[1];
	my $fsa_filename		= $_[0];
	my $i					= 1;
	my $fsa_seq_count		= `grep -c '>' $fsa_filename`;
	my $vector_line_count	= `wc $vector_filename | awk '{print \$1}'`;

	# Test if number of vectores/lines is the same in TAB file as in FASTA file
	unless ($vector_line_count == $fsa_seq_count)	{ 
		printf ("%-40s\t==>\t%-50s\n", "# ERROR checkFileTypeVECTOR", "VECTOR file, $vector_filename, does not have the correct number of lines");				
		$status = "FALSE"; return $status; }

	# Open the VECTOR file and test each vector for problems
	open(FP0,$vector_filename);

	while( my $line = <FP0> )  {
		my @elements = split("\t", $line);
		
#		unless (scalar(@elements) == 82)	{ 
		unless (scalar(@elements) == 79)	{ 
			printf ("%-40s\t==>\t%-50s\n", "# ERROR checkFileTypeVECTOR", "VECTOR, $vector_filename, in line $i does not have the right number of elements");	
			$status = "FALSE"; return $status; }

		foreach my $item (@elements)		{ if ($item eq '') { 
			printf ("%-40s\t==>\t%-50s\n", "# ERROR checkFileTypeVECTOR", "VECTOR, $vector_filename, in line $i has empty elements");							
			$status = "FALSE"; return $status; }}

		unless ($elements[0] =~ /[a-zA-Z]/)	{ 
			printf ("%-40s\t==>\t%-50s\n", "# ERROR checkFileTypeVECTOR", "VECTOR, $vector_filename, element 1 is a number");									
			$status = "FALSE"; return $status; }
		$i++;
	}
	close(FP0);
	return $status;
}

#.............. checkFileTypeTAB ..............
sub checkFileTypeTAB {
	# head -n 1 trial.fsa.tab | grep -o "," | wc
	# head -n 1 trial.fsa.tab | grep -P -o "\t" | wc
	my $status			= "TRUE";
	my $tab_filename	= $_[1];
	my $fsa_filename	= $_[0];
	my $i				= 1;
	my $fsa_seq_count	= `grep -c '>' $fsa_filename`;
	my $tab_line_count	= `wc $tab_filename | awk '{print \$1}'`;

	# Test if number of vectores/lines is the same in TAB file as in FASTA file
	unless ($tab_line_count == $fsa_seq_count+1)	{ 
		printf ("%-30s\t==>\t%-50s\n",  "# ERROR checkFileTypeTAB", "TAB file, $tab_filename, does not have the correct number of lines");					
		$status = "FALSE"; return $status; }

	# Open the TAB file and test each vector for problems
	open(FP0, $tab_filename);
	while( my $line = <FP0> )  {
		my @elements = split("\t", $line);
		next if $line =~ m/number/;

		unless (scalar(@elements) == 88)		{ 
			printf ("%-40s\t==>\t%-50s\n", "# ERROR checkFileTypeTAB", "TAB vector, $tab_filename, in line $i does not have the right number of elements");	
			$status = "FALSE"; return $status; }

		foreach my $item (@elements)			{ if ($item eq '') { 
			printf ("%-40s\t==>\t%-50s\n", "# ERROR checkFileTypeTAB", "TAB vector, $tab_filename, in line $i has empty elements");							
			$status = "FALSE"; return $status; }}

		unless ($elements[2] eq $fsa_filename)	{ 
			printf ("%-40s\t==>\t%-50s\n", "# ERROR checkFileTypeTAB", "TAB vector, $tab_filename, element 3 is not the same as FASTA filename");			
			$status = "FALSE"; return $status; }
		$i++;
	}
	close(FP0);
	return $status;
}

#.............. checkFileTypeFASTA ..............
sub checkFileTypeFASTA {
	# Test the first 10 lines of the file
	my ($filename)	= @_;
	my $i			= 0;
	my @lines;

	open(FP0, $filename);

	while( my $line = <FP0> )  {
		next if ($line =~ /^\s*$/);
		last if ($i > 10);
		chomp $line;	 
		$line =~ s/\s+//g;
		push(@lines, $line);
		$i++;
	}
	close FP0;
	my $filetype	= "nonFASTA";   
	if ( @lines >= 2 && $lines[0] =~ /^>\w/ && $lines[1] =~ /^[a-zA-Z]+$/ ) { $filetype		= "FASTA" }
	return $filetype;
}

#=======================================================================================
#=======================================================================================
#===========================	PATH AND FILE CHECK SUBS	============================
#=======================================================================================
#=======================================================================================

#.............. checkPATHS ..............
sub checkPATHS {
	unless (-d $pybrain_path)	{	#.............. Pybrain path ..............
		printf ("%-30s:\t%-50s\n", "# ERROR checkPATHS", "nothing called $pybrain_path or is not a directory"); 
		printf ("%-30s:\t%-50s\n", "# ERROR checkPATHS", "Path used by CMGfunc_testNetworks.py");
		printf ("%-30s:\t%-50s\n", "# ERROR checkPATHS", "Python pybrain library path: $pybrain_path");	
		exit }

	unless (-d $ProtFun_path)	{	#.............. ProtFun path ..............
		printf ("%-30s:\t%-50s\n", "# ERROR checkPATHS", "nothing called $ProtFun_path or is not a directory"); 
		printf ("%-30s:\t%-50s\n", "# ERROR checkPATHS", "Path used by CMGfunc_calcFeatures.pl");
		printf ("%-30s:\t%-50s\n", "# ERROR checkPATHS", "Analysis_directory : $ProtFun_path");		
		exit }

	unless (-d $CMGfunc_path)	{	#.............. CMGfunc path ..............
		printf ("%-30s:\t%-50s\n", "# ERROR checkPATHS", "nothing called $CMGfunc_path or is not a directory"); 
		printf ("%-30s:\t%-50s\n", "# ERROR checkPATHS", "Path used by CMGfunc.pl");
		printf ("%-30s:\t%-50s\n", "# ERROR checkPATHS", "CMGfunc_directory : $CMGfunc_path");		
		exit }

	unless (-d $NET_path)		{	#.............. CMGfunc path ..............
		printf ("%-30s:\t%-50s\n", "# ERROR checkPATHS", "nothing called $NET_path or is not a directory"); 
		printf ("%-30s:\t%-50s\n", "# ERROR checkPATHS", "Path used by CMGfunc.pl");
		printf ("%-30s:\t%-50s\n", "# ERROR checkPATHS", "NET_directory : $NET_path");				
		exit }
}

#.............. checkSCRIPTS ..............
sub checkSCRIPTS {
	#unless (-e $CMGfunc_path."/CMGfunc_analyzeResults.pl")	{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "nothing called $CMGfunc_path/CMGfunc_analyzeResults.pl");exit }
	unless (-e $CMGfunc_path."/CMGfunc_calcFeatures.pl")	{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "nothing called $CMGfunc_path/CMGfunc_calcFeatures.pl"); 	exit }
	unless (-e $CMGfunc_path."/CMGfunc_normFeatures.pl")	{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "nothing called $CMGfunc_path/CMGfunc_normFeatures.pl"); 	exit }
	unless (-e $CMGfunc_path."/CMGfunc_testNetworks.py")	{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "nothing called $CMGfunc_path/CMGfunc_testNetworks.py"); 	exit }

	unless (-e $ProtFun_path."/ape")	{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no file called $ProtFun_path/ape"); 		exit }
	unless (-e $ProtFun_path."/nmerge")	{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no file called $ProtFun_path/nmerge");	exit }
	unless (-e $ProtFun_path."/nseg")	{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no file called $ProtFun_path/nseg");		exit }
	unless (-e $ProtFun_path."/pseg")	{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no file called $ProtFun_path/pseg");		exit }
	unless (-e $ProtFun_path."/pmerge")	{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no file called $ProtFun_path/pmerge");	exit }
	unless (-e $ProtFun_path."/prop")	{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no file called $ProtFun_path/prop");		exit }

	unless (-e $ProtFun_path."/protparam.py")			{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no file called $ProtFun_path/protparam.py"); 		exit }
	unless (-e $ProtFun_path."/psa2msa")				{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no file called $ProtFun_path/psa2msa"); 				exit }
	unless (-e $ProtFun_path."/runpsipred")				{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no file called $ProtFun_path/runpsipred"); 			exit }
	unless (-e $ProtFun_path."/runpsipred_single")		{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no file called $ProtFun_path/runpsipred_single");	exit }
	unless (-e $ProtFun_path."/seg")					{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no file called $ProtFun_path/seg"); 					exit }
	unless (-e $ProtFun_path."/signalp")				{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no file called $ProtFun_path/signalp"); 				exit }
	#unless (-e $ProtFun_path."/tmhmm-2.0c/bin/tmhmm")	{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no file called $ProtFun_path/tmhmm-2.0c/bin/tmhmm"); exit }

	unless (-d $ProtFun_path."/ape-1.0")				{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no directory called $ProtFun_path/ape-1.0"); 			exit }
	unless (-d $ProtFun_path."/bio-tools-psort-all")	{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no directory called $ProtFun_path/bio-tools-psort-all");	exit }
	unless (-d $ProtFun_path."/libpsortb-1.0")			{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no directory called $ProtFun_path/libpsortb-1.0"); 		exit }
	unless (-d $ProtFun_path."/pftools")				{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no directory called $ProtFun_path/pftools"); 			exit }
	unless (-d $ProtFun_path."/psortb")					{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no directory called $ProtFun_path/psortb"); 				exit }
	unless (-d $ProtFun_path."/prop-1.0c")				{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no directory called $ProtFun_path/prop-1.0c"); 			exit }
	unless (-d $ProtFun_path."/psipred-2.5")			{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no directory called $ProtFun_path/psipred-2.5"); 		exit }
	unless (-d $ProtFun_path."/signalp-4.1")			{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no directory called $ProtFun_path/signalp-4.1"); 		exit }
	#unless (-d $ProtFun_path."/tmhmm-2.0c")				{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "no directory called $ProtFun_path/tmhmm-2.0c"); 			exit }
}

