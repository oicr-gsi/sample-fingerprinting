#!/usr/bin/perl -w

#
# Script for prepering the data and plotting individual genotype fingerprints of .bam files
#

=head2 SYNOPSIS

create_fin.pl script is for parsing genotype and coverage data into a little .fin file that is
used to produce a B<fin>gerprint image for a file.

create_fin.pl --refvcf=dbsnp137.hg19.vcf.gz --genotype=0000.snps.raw.vcf --coverage=tmp/0000.sample_interval_summary --datadir=tmp/ --outdir=data/finfiles/ --basename=0000

=head2 OPTIONS

B<--refvcf=>   ref vcf file                 - a gzipped with annotaed SNPs (for example, data from dbSNP database) - used to extract coordinates of individual SNPs

B<--genotype=> genotype vcf file            - input SNPs in vcf format, to be comapared with reference data and marked accordingly in .fin file

B<--coverage=> sample interval summary file - coverage file produced with BEDtools

B<--datadir=>  data dir                     - that's where all input files are

B<--outdir=>   output dir (optional)        - optional destination for output files, becomes the same as datadir if not set 

B<--basename=> base name                    - this is to be used as the basename for .fin file

=cut

use strict;
use Getopt::Long;
use constant DEBUG=>0;

my(%refs,@bitstrip);

my $USAGE = "create_fin.pl --refvcf=[ref vcf file]
                           --genotype=[genotype vcf file]
                           --coverage=[sample interval summary file]
                           --datadir=[data dir]
                           --outdir=[output dir] (optional)
                           --basename=[base name]\n";

my($refvcf,$vcf,$summary,$datadir,$snpcount,$basename,$outdir);
my $result = GetOptions ('--refvcf=s'   => \$refvcf,
                         '--genotype=s'=>  \$vcf,
                         '--coverage=s' => \$summary,
                         '--datadir=s'  => \$datadir,
                         '--outdir=s'   => \$outdir,
                         '--basename=s' => \$basename);

die $USAGE if (!$refvcf || !$summary || !$datadir || !$vcf);

$outdir ||= $datadir;
my $flags = {zerocov => "N",   # No coverage
             nosnp   => "M"};  # no SNP, base matches the reference

# Read reference SNPs into a hash

my $refhandle = $refvcf=~/\.gz$/ ? "zcat $refvcf | grep -v \"^#\" | " : "grep -v \"^#\" $refvcf |";
open(REFVCF,$refhandle) or die "Couldn't read from reference vcf file [$refvcf]";
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
 $basename = $' if $basename=~m!/!;

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

open(OUT,">$outdir/$basename.fin") or die "Couldn't write to [$outdir/$basename.fin]";

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

