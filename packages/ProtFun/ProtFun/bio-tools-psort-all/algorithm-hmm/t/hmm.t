# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More;
BEGIN { plan tests => 5 };

BEGIN { use_ok('Algorithm::HMM'); }

ok(1); # If we made it this far, we're ok.

#########################

# Insert your test code below, the Test module is use()ed here so read
# its man page ( perldoc Test ) for help writing this test script.

print("Creating new Algorithm::HMM\n");
my $hmm = new Algorithm::HMM(Model => 'sample.hmm');
ok(ref($hmm) ne "", "Ensuring we could load a model");

print("Saving model...\n");
ok(($hmm->save("sample.hmm.1") == 0), "Testing saving a model");

print("Loading saved model...\n");
ok(($hmm->load("sample.hmm.1") == 0), "Testing reloading a model");
