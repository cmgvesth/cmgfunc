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
Pipeline for constructing artificial neural network from FASTA file of protein sequences

#--------------------------------------------------------------
# Input source
#--------------------------------------------------------------

#--------------------------------------------------------------
# Input
#--------------------------------------------------------------
proteins.fsa : FASTA file of protein sequences
head description.fsa
>9832486_PF00004_23
VTSPQSAATDARNAMVAGMLASGISVNGLQPSHNPQVALQMFTTATTIDPSMCDAWLARVLAGDSDIEVLANAWATVRTFGWETRRLGLSDLEFHPEVSDGMFLRLRITSVESLAAAYAAALAEAKRYEEAAKLLDGIEPRNPFETELVGYVRGMLYFRTGRWPDVLKQFPAAAPWRQPELKAAAAAMATTALASLGVFEEATRRAQEAIEGDRVPGATNVALYTQGMCMRHLGREDEAAELLRRVYSRDAKFAPAREALDNLDYRLILTDPETIEARTDPWDPDSAPTREQAEAHRHAEEAARYLAEGDAELAAMLGMEQAKREIKLIKATTKVNLARTKMGLPVPVTSRHTLLLGPPGTGKTSVARAFTKQLCGLTVLRKPLVVETNRSKLLGRYMADAEKNTEELLEGALGGAVFFDEMHTLHERGYSQGDAYGNAIINTLLLYMENHRDELVVFGAGYANAMDKMLEVNQGLRRRFSTVIEFYSYTPPELLELTKMMGAENEDVVTDEAVEGLLPSYTKFYADENYSEDGDLIRGIDVLGNAGFVRNVVEKARDHRSFRLDDSELDAVLAGDVTEFSDEWLHKFKELTREDLAEGLSAAVAEKKPS
>9892486_PF00004_23
MTRPQAAAEDARNAMVAGLLASGISVNGLQPSHNPQVAAQMFTTATRLDPKMCDAWLARLLAGDQSIEVLAGAWAAVRTFGWETRRLGVTDLQFRPEVSDGLFLRLAITSVDSLACAYAAVLAEAKRYQEAAELLDATDPRHPFDAELVSYVRGVLYFRTKRWPDVLAQFPEATQWRHPELKAAGAAMATTALASLGVFEEAFRRAQEAIEGDRVPGAANIALYTQGMCLRHVGREEEAVELLRRVYSRDAKFTPAREALDNPNFRLILTDPETIEARTDPWDPDSAPTRAQTEAARHAEMAAKYLAEGDAELNAMLGMEQAKKEIKLIKSTTKVNLARAKMGLPVPVTSRHTLLLGPPGTGKTSVARAFTKQLCGLTVLRKPLVVETSRTKLLGRYMADAEKNTEEMLEGALGGAVFFDEMHTLHEKGYSQGDPYGNAIINTLLLYMENHRDELVVFGAGYAKAMEKMLEVNQGLRRRFSTVIEFFSYTPQELIALTQLMGRENEDVITEEESQVLLPSYTKFYMEQSYSEDGDLIRGIDLLGNAGFVRNVVEKARDHRSFRLDDEDLDAVLASDLTEFSEDQLRRFKELTREDLAEGLRAAVAEKKTK

description.txt: Text information:
#Short name	Description
Type of hydrolase	This hydrolase was found to perform very well undercondition X Y and X which is uniq...... 

#--------------------------------------------------------------
# Output
#--------------------------------------------------------------
HTML network file
Text file

#--------------------------------------------------------------
# USAGE
#--------------------------------------------------------------
perl CMGfunc_createNetwork.pl -fasta file.fsa -txt description.txt
=cut

#-------- SET PATHS ------
my $defaultNet_path = "/home/cmgfunc/CMGfunc/NETS/LOCAL";
unless (-d $defaultNet_path) { printf ("%-40s:\t%-50s\n", "# ERROR checkPATHS", "defaultNet_directory : $defaultNet_path");	exit }

my $CMGfunc_path	= "/home/cmgfunc/CMGfunc/indirect";
my $pybrain_path	= "/home/cmgfunc/pybrain";
my $ProtFun_path	= "/home/cmgfunc/ProtFun";
my $NET_path 		= $defaultNet_path;
#-------- GET OPTIONS ------

print "#==================================================================================\n";

#.............. Get options ..............
my ($fasta, $txt);
&GetOptions ("fasta:s" =>  \$fasta, "txt:s" => \$txt);

#.............. Input files and paths exists ..............
if (!defined $fasta)	
{ printf ("%-30s:\t%-50s\n", "# ERROR", "Protein FASTA file not defined, example: perl CMGfunc_createNetworks.pl -fasta description.fsa");	&usage;	exit; }
if (!defined $txt)		
{ printf ("%-30s:\t%-50s\n", "# ERROR", "Text description file not defined, example: perl CMGfunc_createNetwork.pl -txt description.txt");	&usage;	exit; }



#.............. Test that paths exists or exit ..............
&checkPATHS();

#.............. Test that scripts exists in the right path or exit ..............
&checkSCRIPTS();

#.............. Print usage ..............
sub usage {	
	printf ("%-30s:\t%-50s\n", "# USAGE", "perl CMGfunc_createNetworks.pl -fasta <name of FASTA file> -txt <name of description file>");	
	printf ("%-30s:\t%-50s\n", "# EXAMPLE", "perl CMGfunc_createNetworks.pl -fasta description.fsa -txt description.txt");
	printf ("%-30s:\t%-50s\n", "# INFO", "First name of FASTA and text file must be the same as this links the two files together");
	exit( 1 );	}


#.............. Test that file is protein FASTA or exit ..............
unless (&checkFileTypeFASTA($fasta) eq "FASTA")						{	printf ("%-30s:\t%-50s\n", "# ERROR", "File, $fasta, is not FASTA format");	exit }
unless ( ((split("\\.", $fasta))[0]) eq ((split("\\.", $txt))[0]) ) {	printf ("%-30s:\t%-50s\n", "# ERROR", "FASTA and TXT does not have same filename");	exit }

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
printf ("%-30s:\t%-60s|\n", '# Python pybrain library path',$pybrain_path);
printf ("%-30s:\t%-60s|\n", '# Analysis directory path', 	$ProtFun_path);
printf ("%-30s:\t%-60s|\n", '# Path to network files', 		$NET_path);
print "#==================================================================================\n";

#=======================================================================================
#=======================================================================================
#=======================================	MAIN	================================
#=======================================================================================
#=======================================================================================

#-------- CALC. FEATURES ------
unless (&checkFileTypeTAB($fasta, $fasta.".tab") eq "TRUE") 			{ &calcFeatures($fasta) }
unless (&checkFileTypeTAB($fasta, $fasta.".tab") eq "TRUE") 			{ exit }

#-------- NORM. FEATURES ------		
&normFeatures($fasta);
unless (&checkFileTypeVECTOR($fasta, $fasta.".tab.vec") eq "TRUE")		{ exit }

#-------- CREATE TRAIN AND TEST ------
if (-e $fasta.".tab.vec.train") {
	printf ("%-40s\t==>\t%-50s\n", "# ERROR", "$fasta.tab.vec.train already found, delete and re-run");		
	exit;
}
&createTestTrain($fasta) unless (-e $fasta.".tab.vec.train");		

#-------- TRAIN NETWORKS ------
if (glob($defaultNet_path."/".$fasta."*.local.net")) {
	my $file = glob($defaultNet_path."/".$fasta."*.local.net"); 
	printf ("%-40s\t==>\t%-50s\n", "# ERROR", "$file already found, delete and re-run");		
	exit;
}
&trainNetworks($fasta) unless glob($fasta."*.local.net");		

`mkdir $defaultNet_path` unless (-d $defaultNet_path );
foreach (glob($fasta."*.local.net*")) {
	`mv $_ $defaultNet_path`;
}
`mv $txt $defaultNet_path`;
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
	unless (-e $fasta.".tab")	{ printf ("%-40s\t==>\t%-50s\n", "# ERROR", "$CMGfunc_path/CMGfunc_calcFeatures.pl did not return file $fasta.tab");exit }
}

#.............. normFeatures ..............
sub normFeatures {
	printf ("%-40s\t==>\t%-50s\n", "# RUNNING CMGfunc_normFeatures.pl", "will create file: $fasta.tab.vec");
	system ("perl $CMGfunc_path/CMGfunc_normFeatures.pl -tab $fasta.tab");
	unless ($? == 0)				{ printf ("%-40s\t==>\t%-50s\n", "# ERROR", "$CMGfunc_path/CMGfunc_normFeatures.pl did not run with success");				exit }
	unless (-e $fasta.".tab.vec")	{ printf ("%-40s\t==>\t%-50s\n", "# ERROR", "$CMGfunc_path/CMGfunc_normFeatures.pl did not return file $fasta.tab.vec");	exit }
}

#.............. createTestTrain ..............
sub createTestTrain {
	printf ("%-40s\t==>\t%-50s\n", "# RUNNING CMGfunc_createTestTrain.sh", "will create files: $fasta.tab.vec.test and $fasta.tab.vec.train");
	system ("$CMGfunc_path/CMGfunc_prepare_train.sh $fasta.tab.vec $CMGfunc_path/neg.tab");
	unless ($? == 0)					{ printf ("%-40s\t==>\t%-50s\n", "# ERROR", "$CMGfunc_path/CMGfunc_createTestTrain.sh did not run with success");					exit }
	unless (-e $fasta.".tab.vec.train")	{ printf ("%-40s\t==>\t%-50s\n", "# ERROR", "$CMGfunc_path/CMGfunc_createTestTrain.sh did not return file $fasta.tab.vec.train");	exit }
	unless (-e $fasta.".tab.vec.test")	{ printf ("%-40s\t==>\t%-50s\n", "# ERROR", "$CMGfunc_path/CMGfunc_createTestTrain.sh did not return file $fasta.tab.vec.test");	exit }

}

#.............. trainNetworks ..............
sub trainNetworks {
	printf ("%-40s\t==>\t%-50s\n", "# RUNNING CMGfunc_trainNetworks.py", "will create file: $fasta.local.net");
	system ("python $CMGfunc_path/CMGfunc_trainNetworks.py -i $fasta.tab.vec.train -t $fasta.tab.vec.test");
	unless ($? == 0)				{ printf ("%-40s\t==>\t%-50s\n", "# ERROR", "$CMGfunc_path/CMGfunc_trainNetworks.py did not run with success");			exit }
	unless (glob($fasta."*.local.net"))	{ printf ("%-40s\t==>\t%-50s\n", "# ERROR", "$CMGfunc_path/CMGfunc_trainNetworks.py did not return file $fasta.local.net");	exit }
}


#=======================================================================================
#=======================================================================================
#===========================	FILETYPE SUBRUTINES	====================================
#=======================================================================================
#=======================================================================================

#.............. checkFileTypeVECTOR ..............
sub checkFileTypeVECTOR {
	my $status				= "TRUE";
	my $vector_filename 	= $_[1];
	my $fsa_filename		= $_[0];
	my $i					= 1;
	my $fsa_seq_count		= `grep -c '>' $fsa_filename`; chomp $fsa_seq_count; $fsa_seq_count =~ s/\n//;
	my $vector_line_count	= `wc $vector_filename | awk '{print \$1}'`;  chomp $vector_line_count; $vector_line_count =~ s/\n//;

	# Test if number of vectores/lines is the same in TAB file as in FASTA file
	unless ($vector_line_count == $fsa_seq_count)	{ 
		printf ("%-40s\t==>\t%-50s\n", "# ERROR checkFileTypeVECTOR", "VECTOR file, $vector_filename, does not have the correct number of lines");				
		$status = "FALSE"; return $status; }

	# Open the VECTOR file and test each vector for problems
	open(FP0,$vector_filename);

	while( my $line = <FP0> )  {
		my @elements = split("\t", $line);
		#print scalar(@elements),"\n";
		unless (scalar(@elements) == 79)	{ 
			printf ("%-40s\t==>\t%-50s\n", "# ERROR checkFileTypeVECTOR", "VECTOR, $vector_filename, in line $i does not have the right number of elements");	
			$status = "FALSE"; return $status; }

		foreach my $item (@elements)		{ if ($item eq '') { 
			printf ("%-40s\t==>\t%-50s\n", "# ERROR checkFileTypeVECTOR", "VECTOR, $vector_filename, in line $i has empty elements");							
			$status = "FALSE"; return $status; }}

		$i++;
	}
	close(FP0);
	return $status;
}

#.............. checkFileTypeTAB ..............
sub checkFileTypeTAB {
	my $status			= "TRUE";
	my $tab_filename	= $_[1];
	my $fsa_filename	= $_[0];
	my $i				= 1;
	my $fsa_seq_count	= `grep -c '>' $fsa_filename`; chomp $fsa_seq_count; $fsa_seq_count =~ s/\n//;
	my $tab_line_count	= `wc $tab_filename | awk '{print \$1}'`; chomp $tab_line_count; $tab_line_count =~ s/\n//;

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
	unless (-e $CMGfunc_path."/CMGfunc_trainNetworks.py")	{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "nothing called $CMGfunc_path/CMGfunc_trainNetworks.py"); exit }
	unless (-e $CMGfunc_path."/CMGfunc_prepare_train.sh")	{ printf ("%-30s:\t%-50s\n",  "# ERROR checkSCRIPTS", "nothing called $CMGfunc_path/CMGfunc_prepare_train.sh"); exit }

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


