# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More;
BEGIN { plan tests => 2 };

use Bio::Tools::PSort;

ok(1, "PSort module loaded"); # If we made it this far, we're ok.

#########################

# Insert your test code below, the Test module is use()ed here so read
# its man page ( perldoc Test ) for help writing this test script.

print("Creating new Bio::Tools::PSort\n");
my $psort = new Bio::Tools::PSort();
ok(ref($psort) ne "", "PSort object made");
