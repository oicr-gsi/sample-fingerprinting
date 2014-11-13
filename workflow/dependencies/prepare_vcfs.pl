#!/usr/bin/perl -w

# =================================================================================================
# A simple script wraps bgzip and tabix so that we don't have missing apostrophes problem in seqware
# =================================================================================================

use strict;
use Getopt::Long;
my($bgzip,$datadir,$tabix);

my $result = GetOptions ('datadir=s'    => \$datadir, # working (output) directory
                         'bgzip=s'      => \$bgzip,  # directory with temporary GATK files
                         'tabix=s'      => \$tabix);  # file with calculated indexes

opendir(DIR,$datadir) or die "Was unable to read from datadir [$datadir]";
my @files = grep {/vcf$/} readdir(DIR);
closedir DIR;
$datadir.="/" if $datadir!~m!/$!;

map {my $file = $datadir.$_;`$bgzip -c $file > $file.gz && $tabix -p vcf $file.gz`;} @files;
print STDERR "Finished preparing vcf files\n";

