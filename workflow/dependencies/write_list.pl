#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Data::Dumper;
use constant DEBUG=>0;
#
# Simple script that takes comma-delimited list of files and converts commas to new lines.
# Print to standard output ================================================================
#
my($datadir,$segments);
my $result = GetOptions ('datadir=s'     => \$datadir,    # working (output) directory with vcf files
                         'segments=s'    => \$segments);  # one or two segments (for composing a merged list) i.e. 0:49,100:149

my $USAGE = "write_list --datadir [dir with vcf.gz files] --segments [optional, segments of file list to make a merged list]\n";
die $USAGE if (!$datadir);

opendir(DIR,"$datadir") or die "Couldn't read from directory [$datadir]";
my @entries = grep{/.vcf.gz$/} grep {!/.vcf.idx/} sort readdir DIR; 
closedir(DIR);

my @segments;

if ($segments) {
 my @chunks = split(",",$segments);
 my $end = @chunks < 2 ? scalar(@chunks) : 2; # Do not accept more than one segment
 for (my $i=0; $i < $end; $i++) {
  my @coords = split(":",$chunks[$i]);
  if (@coords == 2) {
   push(@segments,[$coords[0],$coords[1]]);
  }
 }
} else {
 @segments = ([0,$#entries]);
}

my @vcfs = map{$_ if(/\.vcf.gz$/)} @entries;
print STDERR "Got ".scalar(@vcfs)." vcfs\n" if DEBUG;
print STDERR Dumper(@segments) if DEBUG;

my $split = 0;
foreach my $s (@segments) {
  VCF:
  for (my $v=0; $v<@vcfs; $v++) {
     if ($vcfs[$v]!~/\S/){next VCF;}
     if ($v >= $s->[0] && $v <= $s->[1]) {
       $vcfs[$v]=~s/\.snps\.raw\.vcf.*//;
       print $vcfs[$v]."\n";
     }
  }
  print "\n" if 0 == $split++;
}

