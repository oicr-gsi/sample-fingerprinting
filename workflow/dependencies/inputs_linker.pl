#!/usr/bin/perl -w
#==================================================================================
#
# Script for linking inputs into one datadir (needed for Sample Fingerprinting 2.0)
#
#==================================================================================

use strict;
use Data::Dumper;
use Getopt::Long;

use constant DEBUG=>1;
my($list,$datadir,$extensions);
my(@extensions,@basepaths);
my $result = GetOptions ('list=s'            => \$list, #Comma-separated list with basenames (full path sans extension), extension should be passed in 'extensions' 
                         'datadir=s'         => \$datadir, # directory for symbolic links to inputs
                         'extensions=s'      => \$extensions); #Comma-separated extensions of all files that need to be linked as well

# split extensions
@extensions = split(",",$extensions);
@basepaths  = split(",",$list);

# for each path, for each extension link
for (my $f = 0; $f < @basepaths; $f++) {
    for (my $e = 0; $e < @extensions; $e++) {
      my $file = $basepaths[$f].$extensions[$e];
      print STDERR "Linking $file into $datadir\n" if DEBUG;
      `ln -s $file -t $datadir`;
    }
}

