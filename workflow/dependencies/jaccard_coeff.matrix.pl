#!/usr/bin/perl -w

# =====================================================
# calculate jaccard cefficients for a list of vcf files
# output as giant matrix (Version for workflow)
# =====================================================

use Getopt::Long;
use Data::Dumper;
use strict;

use constant DEBUG =>1;
# Below is the dafault for vcf_compare, should not be used when workflow runs
my $vcf_compare = "vcftools/bin/vcf-compare";

my(%ids,%files,%matrix,%seen,$list,$studyname,$vcf_path,$oldmatrix,$datadir);
my(%snps,%cols,%scols); # For assigning colors and number of snps (scols for assigning a color to a sample (like PCSI_006)
my @colors = qw/red orange yellow green lightblue blue purple darkgreen brown black/;

my $USAGE="jaccard_coef.matrix.pl --list=[req] --studyname=[req] --existing_matrix=[optional] --vcf-compare=[] ";

my $result = GetOptions ('list=s'            => \$list, # list with filenames
                         'existing_matrix=s' => \$oldmatrix, # file with previously calculater indexes
                         'vcf_compare=s'     => \$vcf_path, # path to vcftools vcf-compare script
                         'studyname=s'       => \$studyname);

if (!$list || !$studyname) {die $USAGE;}
$vcf_compare = $vcf_path if $vcf_path;

# ==========================================================
# First, if we have exisiting matrix, read values from there
# ==========================================================

my @old_matrices = $oldmatrix ? grep{/\S+/} split(",",$oldmatrix) : undef;
if (@old_matrices && @old_matrices > 0) {
foreach my $old (@old_matrices) {
 if (-e $old) {
  print STDERR "We have some older data, will append our results to the existing matrix\n" if DEBUG;
  open(OLD,"<$old") or die "Couldn't read from the file with previously calculated indexes";
  my $idstring = <OLD>;
  chomp($idstring);

  my @tempids = grep {/\D+/}  split("\t",$idstring);
  my $has_snps = $tempids[$#tempids] =~ /SNP/ ? 1 : 0;
  pop @tempids if $tempids[$#tempids] =~ /SNP/;

  map {&id_file($_)} (@tempids);

  while (<OLD>) {
   chomp;
   my @temp = split("\t");
   $temp[0] = $` if $temp[0] =~/.snps.raw.vcf.gz$/;

   for (my $t = 0;$t < @tempids;$t++) {
     $tempids[$t] = $` if $tempids[$t]=~/.snps.raw.vcf.gz$/;
     $seen{$files{$temp[0]}}->{$files{$tempids[$t]}}++;
     $seen{$files{$tempids[$t]}}->{$files{$temp[0]}}++;

     $matrix{$files{$tempids[$t]}}->{$files{$temp[0]}} = $matrix{$files{$temp[0]}}->{$files{$tempids[$t]}} = $temp[$t+1] if ($files{$tempids[$t]} && $files{$temp[0]});
    }

   $snps{$files{$temp[0]}} = $has_snps ? $temp[$#temp] : &calculate_snps($temp[0]);
  }
  close OLD;
 }
}
print STDERR "Identified ".scalar(keys %seen)." interactions and number of SNPs for ".scalar(keys %snps)." files using old matrices\n" if DEBUG;
}


# =================================================================================
# Collect information on files, build matrix using new (and old, if available) data
# =================================================================================

if ($list=~/\,/) {
 # We have a comma-delimited list of files
 my @files = split(",",$list);
 map{&id_file($_)} @files; 
} else {
 open(LIST,"<$list") or die "Cannot read from the list file [$list]";
 map{chomp;&id_file($_)} (<LIST>);
 close LIST;
}

print STDERR scalar(keys %ids)." genotypes collected\n" if DEBUG;

my $count = 1;
foreach my $id(keys %ids) {
 # Take care of SNP number calculation
 $snps{$id} ||= &calculate_snps($ids{$id});
  
 my $sample = $ids{$id};

 print STDERR "Working on ".$count++." of ".scalar(keys %ids)." samples\n";
 SM:
 foreach my $s(keys %ids) {
    if ($s eq $id) {
   $matrix{$id}->{$s} = 1;
   next SM;
  }
  
  my $file1 = $ids{$id};
  my $file2 = $ids{$s};

  if ($seen{$id}->{$s}) {next;} 
  $seen{$id}->{$s}++;
  $seen{$s}->{$id}++;
  
  $file1 .= ".snps.raw.vcf.gz" if $file1 !~/gz$/;
  $file2 .= ".snps.raw.vcf.gz" if $file2 !~/gz$/;
  print STDERR "Will check files $file1 and $file2\n" if DEBUG;
  if (! -e $file1 || ! -e $file2){print STDERR "File(s) $file1 or $file2 not FOUND!";next;}
  print STDERR "Will run $vcf_compare $file1 $file2 ...\n" if DEBUG;
  my @compares = `$vcf_compare $file1 $file2 | grep \"^VN\"`;
  my @numbers = (0,0,0);

  for (my $i=0; $i<@compares; $i++) {
    if ($compares[$i] =~ /$ids{$id}/ && $compares[$i] =~ /$ids{$s}/) {
     $numbers[2] = $1 if ($compares[$i]=~/\t(\d+)\t/);
    } elsif ($compares[$i] =~ /$sample/) {
     $numbers[0] = $1 if ($compares[$i]=~/\t(\d+)\t/);
    } elsif ($compares[$i] =~ /$ids{$s}/) {
     $numbers[1] = $1 if ($compares[$i]=~/\t(\d+)\t/);
    }
  }
 
  map{chomp} @compares;
  my $union = 0;
  map{$union+=$numbers[$_]} (0..2);

  $matrix{$id}->{$s} = $union > 0 ? sprintf "%.3f",$numbers[2]/$union : print 0; 
 }
}

# ====================================
# Printing out scores in matrix-style
# ====================================

my @heads;
map {/.snps.raw.vcf.gz/ ? push(@heads,$`) : push(@heads,$ids{$_})} (sort keys %matrix);
print join("\t",("",@heads,"SNPs")); #,"Color","SNPs"));
print "\n";

 foreach my $sample(sort keys %matrix) {
 print $sample=~/.snps.raw.vcf.gz/ ? $` : $ids{$sample};
 TF:
 foreach my $ss(sort keys %matrix) {
   my $value = $matrix{$sample}->{$ss} || $matrix{$ss}->{$sample};
   print $value ? "\t$value" : "\t0";
 }
 # Colors will be assigned using another script
 print $snps{$sample} ? "\t$snps{$sample}" : "\t0";
 print "\n";
 }

# ==================================================
# Subroutine for processing (registering) a vcf file
# ==================================================

sub id_file {
 my $file = shift @_;
 return if $file!~/\D+/; 
 $file = $` if $file =~/.snps.raw.vcf.gz$/;

 if ($file=~/(\d+)_($studyname.\d+_)/) {
  #print STDERR "Studyname $studyname detected\n" if DEBUG;
  my $id = $2.$1;
  $ids{$id} = $file;
  $files{$file} = $id;
 } elsif ($file=~/($studyname.{0,1}\d+)\.(\S+)/) {
  print STDERR "Studyname $studyname detected\n" if DEBUG;
  my $id = join("_",($1,$2));
  $ids{$id} = $file;
  $files{$file} = $id;
 } else {
  #print STDERR "Studyname $studyname NOT detected\n" if DEBUG;
  $ids{$file} = $file;
  $files{$file} = $file;
 }
}

# ==================================================
# Subroutine for calculating number of SNPs
# ==================================================

sub calculate_snps {
 my $file = shift @_;
 print STDERR "Calculating SNPs for $file\n" if DEBUG;
 if ($file !~/.snps.raw.vcf.gz$/) {
   $file.=".snps.raw.vcf.gz";
 }

 my $result = `zcat $file | grep -v ^# | wc -l`;
 chomp($result);
 return $result;
}
