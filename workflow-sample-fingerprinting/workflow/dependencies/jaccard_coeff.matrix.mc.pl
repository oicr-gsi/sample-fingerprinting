#!/usr/bin/perl -w

=head2 jaccard_coeffi.matrix.mc.pl
  
 aka matrix_maker.pl
 a script that operated on fin files and creates a jaccard coefficient matrix
 
 ./jaccard_coeffi.matrix.mc.pl --list [list of .fin file] --out [output file] --use-uncovered [Optional flag to include locations with N flag]

=cut

use strict;
use Getopt::Long;
use constant DEBUG=>0;

my($list,$findir,$outfile,$use_uncovered,$short);
my $USAGE = "./jaccard_coeffi.matrix.mc.pl --list [list of .fin file] --use-uncovered [Optional flag to include locations with N flag]";
my $results = GetOptions("list=s"        => \$list,
                         "findir=s"      => \$findir,
                         "short-name"    => \$short, # Shorten ids (false by default)
                         "use-uncovered" => \$use_uncovered);

if ((!$list && !$findir)) {die $USAGE;}
=head2 Assembling a list of files
 We read from either list or
 directory containing .fin files
=cut

my @finfiles = ();
if ($list && -e $list) {
  open(FINS,"<$list") or die "Couldn't read from list of finfiles";
  map {chomp;push @finfiles, $_} (<FINS>);
  close FINS;
} elsif (-d $findir) {
  opendir(DIR,"$findir") or die "Couldn't read from directory [$findir]";
  my @fins = grep {/.fin$/} sort readdir DIR;
  map{push @finfiles,join("/",($findir,$_))} @fins;
  closedir DIR;
} else { die "Wrong source of fin files, won't run!";}

my %snps   = ();
my %counts = ();

=head2
 Reading from list of files
 Here we extract the file paths from a list, check if the file exists
 and then use the values to populate snps hash
=cut

foreach my $file(@finfiles) {
 next if $file=~/^#/;
 chomp($file);
 if ( !-e $file || -S $file) {
   print STDERR "File[$file] is absent or empty, won't use\n";
   next;
 }
 my @lines = `tail -n +2 $file | cut -f 3,5`;
 my $id = $file;
 $id=~s!.*/!!;
 if ($short) {
  $id=~s!SWID_\d+_!!;
  $id=~s/(_TS)_.*/$1/;
 } else {
  $id=~s/.fin$//;
 }

 map{chomp;my @temp=split("\t");$snps{$id}->{$temp[0]} = $temp[1];$counts{$id}++ if $temp[1]=~/[ACTG]/} @lines;
}

=head2 Calculate coefficients
 Making Jaccard coefficient matrix
 Depending on use_uncovered flag we either use or don't use SNPs 
 where coverage in one of the samples is missing

 AGTC = Change (SNP called) M = match with the reference N = No coverage

=cut

my %matrix = ();
my $count = 1; #Progress counter
foreach my $s (keys %snps) {
   next if !$counts{$s};
   print STDERR $count++." of ".scalar(keys %snps);
   print STDERR "\r";
   
  foreach my $ss (keys %snps) {
    next if !$counts{$ss};
    if ($s eq $ss) {$matrix{$s}->{$ss} = 1;
                    next;
    }
  next if $matrix{$s}->{$ss}; # Avoid analyzing one pair twice
  $matrix{$s}->{$ss} = &jaccard_coeff($snps{$s}, $snps{$ss});
  $matrix{$ss}->{$s} ||= $matrix{$s}->{$ss};
  }
}

=head2 Make matrix
 Printing the results
 This part just outputs the SNP data
=cut

print "\t".join("\t",(sort keys(%matrix)))."\tSNPs\n";
foreach my $s (sort keys %matrix) {
 print $s;
 map{print $matrix{$s}->{$_} ? "\t".$matrix{$s}->{$_} : "\t0"} (sort keys %matrix);
 print "\t".$counts{$s};
 print "\n";
}

=head2 Coefficient function
 Depending on use_uncovered flag we calculate Jaccard coefficient
 as Intersect/Union
=cut

sub jaccard_coeff {
 my ($sample_one, $sample_two) = @_;
 my ($intersect,$union) = (0,0);

 foreach my $rs (keys %{$sample_one}) {
   print STDERR $rs."\n" if DEBUG;
   if (($sample_one->{$rs} eq 'N' || $sample_two->{$rs} eq 'N') && !$use_uncovered) {
       print STDERR "Flag set to ignore uncovered and we have ".$sample_one->{$rs}." and ".$sample_two->{$rs}."\n" if DEBUG;
       next;
   }
   next if $sample_one->{$rs} eq 'N' && $sample_two->{$rs} eq 'N'; # Just skip these, we don't have any info here;
   $union++;
   if ($sample_one->{$rs} eq $sample_two->{$rs}) {
     print STDERR "Have a match\n" if DEBUG;
     $intersect++;
   }
 }
 $union ||=1; # Avoid devision by zero
 return sprintf "%.4f", $intersect/$union;
}
