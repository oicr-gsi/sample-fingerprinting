#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use Test::More 'no_plan';

#
# A [very] simple test script that checks that R is installed and 
#


# Pre-determined value gets returned from the test script
my $result = `Rscript $Bin/rscripts/test_r.r`;

if (!$result || $result ne "OK") {
 ok( 0, "No R installed, the package won't work" );
} else {
 ok( 1, "R present, all is ok" );
}

