# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More;
BEGIN { plan tests => 4 };

BEGIN { use_ok('Bio::Tools::PSort::ModHMM'); }

ok(1, "Loading ModHMM"); # If we made it this far, we're ok.

#########################

# Insert your test code below, the Test module is use()ed here so read
# its man page ( perldoc Test ) for help writing this test script.

my $seq = 'MMSSPHPMSSSRNTPLGVFYSLLACFYWGMVFVIPSMLGNFADLDIVLTRYSVFGICSLITILYKRSNIFKTVPFFLWKKGILWAFLINIAYYFGIAQAVRYSGSAVTVIIAGLAPIAILFYSNIKKKMLSYSFLLSMSGIIVVGIILSNVSEFQSESSSSLPLYLLGLGCVTAATSIWAGYIICNHDFLEQHSEISPDTWCHMLGISSLIICLPLIILGDTFGITHVTRNFLFHTPLSERCLFIVLCSAMGIFSSSRAIAAWNKASLHLSTALLGALLIFEPIFGWILSYLCKREMPSFQEGLGFFLMLGASLCLLLAQKKASEQETPSETLITTESDANJ';
my $hmmfile = "./HMMS/S_TMHMM_0.92b.hmg";
my $repfile = "./HMMS/replacement_letter_multi.rpl";
my $path = "./HMMS/";

print("Running sequence through Bio::Tools::PSort::ModHMM\n");
my $ret = findHelices($seq, $hmmfile, $repfile, $path);
ok(($ret eq 10), "Running sequence through ModHMM");

print("Retreiving TM helix sites\n");
$ret = getHelices($seq, $hmmfile, $repfile, $path);
ok(($ret->count eq 10), "Retreiving TM helix sites");
