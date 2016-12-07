#!/usr/bin/perl

# Downloads requested references from parameters passed on the command line
# Script originally written for the PCAWG BWA MEM workflow: https://github.com/ICGC-TCGA-PanCancer/Seqware-BWA-Workflow
# Modified for use in other workflows

use strict;
use warnings;
use Getopt::Long;

my $link_dir   = "/data";
my $resources="http://s3.amazonaws.com/pan-cancer-data/pan-cancer-reference/genome.fa.gz,http://s3.amazonaws.com/pan-cancer-data/pan-cancer-reference/genome.fa.gz.fai";
my $verbose;
GetOptions ("resources=s" => \$resources,    # numeric
            "link_dir=s"   => \$link_dir,      # string
            "verbose"  => \$verbose)   # flag
  or die("Error in command line arguments\n");
check_tools();
my @words = split /,/, $resources;

if ($verbose) { print "Making directory $link_dir"; }
system("mkdir -p $link_dir");

foreach (@words) {
    download($link_dir,$_)
}

sub download {
  my ($dir, $url) = @_;
  system("mkdir -p $dir");
  $url =~ /\/([^\/]+)$/;
  my $file = $1;

  print "\nDOWNLOADING HEADER:\n\n";
  my $r = `curl -I $url | grep Content-Length`;
  $r =~ /Content-Length: (\d+)/;
  my $size = $1;
  print "\n+REMOTE FILE SIZE: $size\n";
  my $fsize = -s "$dir/$file";
  print "+LOCAL FILE SIZE: $size\n";

  if (!-e "$dir/$file" || -s "$dir/$file" == 0 || -s "$dir/$file" != $size) {
    my $cmd = "wget --no-verbose -c -O $dir/$file $url"; 
    print "\nDOWNLOADING: $cmd\nFILE: $file\n\n";
    my $r = system($cmd);
    if ($r) {
      print "+DOWNLOAD FAILED!\n";
      $cmd = "lwp-download $url $dir/$file";
      print "\nDOWNLOADING AGAIN: $cmd\nFILE: $file\n\n";
      my $r = system($cmd);
      if ($r) {
        die ("+SECOND DOWNLOAD FAILED! GIVING UP!\n");
      }
    }
  } else {
			print "\n$file FOUND LOCALLY - SKIPPING DOWNLOAD\n"
	}
}

sub check_tools {
  if (system("which curl") || (system("which lwp-download") && system("which wget"))) {
    die "+TOOLS NOT FOUND: Can't find curl and/or one of lwp-download and wget, please make sure these are installed and in your path!\n";
  }
}
