#!/usr/bin/perl -w

#
# Script for prepering the data and plotting individual genotype fingerprints of .bam files
#

use strict;
use Getopt::Long;
use constant DEBUG=>0;

my(%refs,@bitstrip);

my $USAGE = "create_fin.pl --refvcf=[ref vcf file]
                           --genotype=[genotype vcf file]
                           --coverage=[sample interval summary file]
                           --datadir=[data dir]
                           --basename=[base name]\n";

my($refvcf,$vcf,$summary,$ccols,$datadir,$pngname,$snpcount,$basename,$workdir);
my $result = GetOptions ('--refvcf=s'   => \$refvcf,
                         '--genotype=s'=>  \$vcf,
                         '--coverage=s' => \$summary,
                         '--datadir=s'  => \$datadir,
                         '--basename=s' => \$basename);

die $USAGE if (!$refvcf || !$summary || !$datadir || !$vcf);

my $flags = {zerocov => "N",   # No coverage
             nosnp   => "M"};  # no SNP, base matches the reference

# Read reference SNPs into a hash

open(REFVCF,"grep -v \"^#\" $refvcf |") or die "Couldn't read from reference vcf file [$refvcf]";
while(<REFVCF>) {
 $snpcount++;
 my @temp = split("\t");
 $refs{$temp[0]}->{$temp[1]} = {flag =>$flags->{zerocov},
                                dbsnp=>$temp[2],
                                ref  =>$temp[3]};
}
print STDERR scalar(keys %refs)." reference points retrieved\n" if DEBUG;
close REFVCF;

# Read from summary file, compose a .fin file

open(SUMMARY,"<$summary") or die "Couldn't read from summary file [$summary]";
 my $count = 0;
 $basename ||= $summary;

 if ($basename =~m!\.sample_interval_summary!) {
  $basename = $`;
 }

print STDERR "Basename is $basename\n" if DEBUG;

while (<SUMMARY>) {
 next if (/^Target/); # skip the first line
 my @temp = split("\t");
 my @coords = split(":",$temp[0]);  
 next if @coords !=2;
 $count++;
 if (($refs{$coords[0]}->{$coords[1]} && $temp[1] > 0)) {
   $refs{$coords[0]}->{$coords[1]}->{flag} =  $flags->{nosnp}; 
 } else {
   $refs{$coords[0]}->{$coords[1]}->{flag} = $flags->{zerocov};
 }
}
print STDERR $count." coordinates got coverage info\n" if DEBUG;
close SUMMARY;

open(OUT,">$datadir/$basename.fin") or die "Couldn't write to [$datadir/$basename.fin]";

my $vcfhandle = $vcf=~/\.gz$/ ? "zcat $vcf | grep -v \"^#\" | " : "grep -v \"^#\" $vcf |";
open(VCF,$vcfhandle) or die "Couldn't read from vcf file [$vcf]";
while (<VCF>) {
 my @temp = split("\t");
 next if (!$refs{$temp[0]}->{$temp[1]});
 if ($temp[4] =~/[ACGT]{1}/i && $temp[4] ne $temp[3]) {
   $refs{$temp[0]}->{$temp[1]}->{flag} = uc($temp[4]);
 } else { $refs{$temp[0]}->{$temp[1]}->{flag} = $flags->{nosnp};}
}

close(VCF);
print OUT join("\t",("CHROM","POS","ID","SNP","FLAG")),"\n";
foreach my $chr (sort keys %refs) {
 foreach my $coord (sort {$a<=>$b} keys %{$refs{$chr}}) {
  my $snp = $refs{$chr}->{$coord}->{flag} eq 'M' ? $refs{$chr}->{$coord}->{ref}.$refs{$chr}->{$coord}->{ref} : $refs{$chr}->{$coord}->{ref}.$refs{$chr}->{$coord}->{flag};
  $snp = "" if $refs{$chr}->{$coord}->{flag} eq 'N';
  
  print OUT join("\t",($chr,$coord,$refs{$chr}->{$coord}->{dbsnp},$snp,$refs{$chr}->{$coord}->{flag})),"\n";
 }
 print STDERR scalar(keys %{$refs{$chr}})." datapoints will be printed out\n" if DEBUG;
 
}

close OUT;

