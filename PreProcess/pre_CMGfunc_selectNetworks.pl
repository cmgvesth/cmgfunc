#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::SeqIO;
use List::Util qw( min max );
use File::Copy;

=put
#=====================================================================
Script takes the MCC result files highest scoring network.
Then the corresponding network files are copied to a NET directory.
The network files, themselfs do not hold the MSE value, as the value is calculated when the networks are TESTED on a TEST file.
It is therefor nessesary to use the performance files

The script assumes that for arch clusters, only set C has been compiled
It also assumes that for clan sets, A, B and C sets have been created

The script works on three directories each containing a specific type of files
HACK NOTE: If the performance file has not run until end and has not returned a MCC line, the mse is evaluated

perl /home/cmgfunc/MySQL/CMGfunc/pre_CMGfunc_selectNetworks.pl -perdir /home/cmgfunc/MySQL/CMGfunc/DATA/SetCLANS/PERFORMANCE -prodir /home/cmgfunc/MySQL/CMGfunc/DATA/SetCLANS/PROGRESS -netdir /home/cmgfunc/MySQL/CMGfunc/DATA/SetCLANS/NETS

perl /home/cmgfunc/MySQL/CMGfunc/pre_CMGfunc_selectNetworks.pl -perdir /home/cmgfunc/MySQL/CMGfunc/DATA/SetCLANS/PERFORMANCE -prodir /home/cmgfunc/MySQL/CMGfunc/DATA/SetCLANS/PROGRESS -netdir /home/cmgfunc/MySQL/CMGfunc/DATA/SetCLANS/NETS

perl /home/cmgfunc/MySQL/CMGfunc/pre_CMGfunc_selectNetworks.pl -perdir /home/cmgfunc/MySQL/CMGfunc/DATA/SetArchClans/PERFORMANCE -prodir /home/cmgfunc/MySQL/CMGfunc/DATA/SetArchClans/PROGRESS -netdir /home/cmgfunc/MySQL/CMGfunc/DATA/SetArchClans/NETS
#=====================================================================
=cut
my ($per_dir, $pro_dir, $net_dir,%dict, $mcc, %root_names);

#.............. Get options ..............
&GetOptions ("perdir:s" =>  \$per_dir, "prodir:s" => \$pro_dir, "netdir:s" => \$net_dir,);


#.............. Input files and paths exists ..............
unless (defined $per_dir and defined $pro_dir and defined $net_dir )		{	
	printf ("%-30s:\t%-50s\n", "# ERROR", "Path to networfiles not defined");			
	printf ("%-30s:\t%-50s\n", "# EXAMPLE", "perl /home/cmgfunc/MySQL/CMGfunc/pre_CMGfunc_selectNetworks.pl -perdir /home/cmgfunc/MySQL/CMGfunc/DATA/SetCLANS/PERFORMANCE -prodir /home/cmgfunc/MySQL/CMGfunc/DATA/SetCLANS/PROGRESS -netdir /home/cmgfunc/MySQL/CMGfunc/DATA/SetCLANS/NETS");	exit; }

#............. List network files ..................
opendir (NETDIR, $net_dir) or die $!;
while (my $netfile = readdir(NETDIR)) {
	chomp $netfile;
	next if ($netfile =~ m/^\./);
	$root_names{(split("\\.", $netfile))[0]} = 1;
}
closedir(NETDIR);

foreach my $r (keys %root_names) {
	chomp $r;
	my ($mcc_A, $mcc_B, $mcc_C, $max, $max_lab) = (0, 0, 0, 0, "");

	if ($net_dir !~ m/Arch/) {

		$mcc_A = `grep MSE $pro_dir/$r.tab.A*progress | grep test | cut -f 6`; chomp $mcc_A;
		$mcc_B = `grep MSE $pro_dir/$r.tab.B*progress | grep test | cut -f 6`; chomp $mcc_B;
		$mcc_C = `grep MSE $pro_dir/$r.tab.C*progress | grep test | cut -f 6`; chomp $mcc_C;

		if ($mcc_A > $mcc_B) { $max = $mcc_A; $max_lab = "A"}
		else { $max = $mcc_B; $max_lab = "B"}
		if ($max < $mcc_C) { $max = $mcc_C, ; $max_lab = "C"}
	}

	else {
		$mcc_C = `grep MSE $pro_dir/$r.tab.C*progress | grep test | cut -f 6`; chomp $mcc_C; $mcc_C =~ s/\n//;
		#print "MCC: $mcc_C\n";
		next if ($mcc_C eq '');
		#if ($mcc_C eq '') { 
		#Score:  0.000973399448694  loop:  193
		#$mse_C = `tail -n 1 $pro_dir/$r.tab.C*progress | awk '{print $2}'`; chomp $mcc_C; $mcc_C =~ s/\n//;
		#}
		$max_lab = "C";
		#print "MCC: $mcc_C\n";
		$max = $mcc_C;
	}

	my $pro_file = `ls $pro_dir/$r.tab.$max_lab*progress`; chomp $pro_file; $pro_file =~ s/\n//;
	my $per_file = `ls $per_dir/$r.tab.$max_lab*performance`; chomp $per_file; $per_file =~ s/\n//;
	my $net_file = `ls $net_dir/$r.tab.$max_lab*net`; chomp $net_file; $net_file =~ s/\n//;

	if ($max < 0.85) {
		copy($pro_file, $pro_file.".lowmax");
		copy($per_file, $per_file.".lowmax");
		copy($net_file, $net_file.".lowmax");
	} else { 
		copy($pro_file, $pro_file.".max");
		copy($per_file, $per_file.".max");
		copy($net_file, $net_file.".max");

	}
}
