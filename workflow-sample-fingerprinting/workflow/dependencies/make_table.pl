#!/usr/bin/perl -w

=head2 make_table.pl
 
 This script relies on data produced by SampleFingerprinting workflow
 It aims to extract data for lane-level bams and output everything for a follow-up analysis with R 
 (and possibly other means). The table that this script produces has the following columns:

 | FLAG | Sample | FlowCell | Lane | Type | SNPs | Score | SimilarTo |
 
 The information comes from two pieces - index.html and jaccard.matrix.csv

=cut

use strict;
use Getopt::Long;
use Data::Dumper;
use IO::File;
use constant DEBUG=>0;
use constant SEP=>"\t";     # Separator
use constant THRESHOLD=>10; # Require that many SNPs to be considered as SimilarTo sample

my($matrix,$index,$outfile,$snp_index,$single_ok);
my %samples = (); # Sample counter
my $USAGE = "make_table.pl --matrix [jaccard_coefficient.matrix.csv] --index [path to index.html] --out [Optional output file for writing table] --singletons [Allow singletons]\n";
my $opts  = GetOptions('matrix=s'     => \$matrix,
                       'index=s'      => \$index,
                       'singletons'   => \$single_ok,
                       'output|out=s' => \$outfile);


my %data = (); # Put all data here, this is keyed by id found in index.html and sim. matrix

=head2 Reading flagged bams
 We open index.html file and if there are flagged bams
 we put them into a hash %flagged
=cut

my $id_string = `grep -i "Flagged as potential" $index | sed 's!/table.*!!'`;
my @temp = split('<tr align="LEFT" valign="BOTTOM"><td>',$id_string);
my @flagged = ();
map{chomp;s/<.*//;if (/\S+/){/_(SWID_\d+)$/ ? push @flagged, join("_",($1,$`)) : push @flagged,$_;}} @temp;

print STDERR scalar(@flagged)." Flagged samples found\n" if DEBUG;
my %flagged = map{$_=>1} @flagged;

map{print STDERR $_."\n"} (@flagged) if DEBUG;

=head2 Reading similarity data
 We open jaccard coefficient matrix this time,
 reading 1) coefficients 2) number of SNPs
 
 we also count how many bams each donor has here
=cut

my $firstline = `head -n 1 $matrix`;
$firstline=~s/^0*//; # Funny entries should be removed by this command
chomp($firstline);

# Get ids as an array
my @heads = split("\t",$firstline);
map {if ($heads[$_]=~/^SNP/){$snp_index = $_}} (0..$#heads);

# Load number of SNPs in a hash
my $snp_column = $snp_index + 1;
my @snp_lines = `cut -f 1,$snp_column $matrix`;
my %snps = map{chomp;my @tmp = split("\t");$tmp[0]=>$tmp[1]} @snp_lines;

# Read from matrix file
my $matrix_pipe = "tail -n +2 $matrix |";
open(MATRIX,$matrix_pipe) or die "Couldn't open pipe for matrix file [$matrix]";
while(my $line = <MATRIX>) {
 chomp($line);
 my @fields = split("\t",$line);
 # This regex assumes that id complies with OICR naming ideosyncrasis
 # | FLAG | Sample | FlowCell | Lane | Type | SNPs | Score | SimilarTo |
 if ($fields[0] =~/^SWID_\d+_(\S+?_\d+)_\S+?_(\S+?)_\S+?_\d+_[A-Z]+?(_\d+)*_(\d+_\S+?_\d+_\S+?)_(\S+?)_L(\d+)_/) {
   $data{$fields[0]} = {FLAG     => $flagged{$fields[0]} ? 1 : 0,
                        Sample   => $1,
                        FlowCell => $4,
                        Lane     => $6,
                        Barcode  => $5,
                        Type     => $2};
   $samples{$1}++;
                        
   my ($score,$simto);
   $score = 0;
   $simto = "NA";
   FIELD:
   for (my $f = 1; $f < $snp_index; $f++) {
     if ($heads[$f] eq $fields[0]) {next FIELD;}
     if ($score < $fields[$f] && $snps{$heads[$f]} >= THRESHOLD) {
         # If sample was not flagged, do not try to match it with unrelated samples
         #if (!$flagged{$fields[0]} && $heads[$f] !~/$data{$fields[0]}->{Sample}/) {next FIELD;} 
         $score = $fields[$f];
         $simto = $heads[$f];
     }
   }
   $data{$fields[0]}->{Score}     = $score;
   $data{$fields[0]}->{SimilarTo} = $simto;
   $data{$fields[0]}->{SNPs}      = $snps{$fields[0]};
 }
}
close(MATRIX);

=head2 Printing
 If needed, we open a file handle
 Otherwise, just print to STDOUT

 Note that we print out singletons only if requested
=cut

my @toprint = qw(FLAG Sample FlowCell Lane Barcode Type SNPs Score SimilarTo);
my $fh = new IO::File();
if ($outfile) {
 $fh->open(">$outfile") or die "Couldn't write to [$outfile]";
 print $fh join("\t",@toprint)."\n";
}else {
 print join("\t",@toprint)."\n";
}

foreach my $id (sort keys %data) {
 my @values = map{$data{$id}->{$_}} @toprint;
 if (scalar(@values) == scalar(@toprint)) {
     if (!$single_ok && $samples{$data{$id}->{Sample}} < 2) {next;} # Filter out singletons if there's no ok flag
     $outfile ? print $fh join(SEP,@values)."\n" : print join(SEP,@values)."\n";
 } else {
   print STDERR "id [$id] misses values!!!\n";
 }
}

if ($outfile) { $fh->close; }

