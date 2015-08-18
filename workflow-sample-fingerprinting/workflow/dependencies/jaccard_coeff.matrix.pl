#!/usr/bin/perl -w

# =====================================================
# calculate jaccard cefficients for a list of vcf files
# output as giant matrix (Version for workflow)
# =====================================================

use Getopt::Long;
use Data::Dumper;
use strict;


use constant DEBUG =>0;
# Below is the dafault for vcf_compare, should not be used when workflow runs
my $vcf_compare = "vcftools/bin/vcf-compare";
my(%ids,@sublists,%files,%matrix,%seen,$list,$studyname,$vcf_path,$oldmatrix,$datadir,$path_to_tabix);
my(%snps); # For assigning number of snps 
my $USAGE="jaccard_coef.matrix.pl --list=[req] --studyname=[req] --datadir=[optional] --existing_matrix=[optional] --vcf-compare=[]";
my $result = GetOptions ('list=s'            => \$list, # list with filenames (empty line may devide subset 1 from subset 2, so sub1 compared to sub2 but not to itself) 
                         'existing_matrix=s' => \$oldmatrix, # file with previously calculater indexes
                         'datadir=s'         => \$datadir, # directory with vcf files
                         'tabix=s'           => \$path_to_tabix, # path to tabix dir
                         'vcf-compare=s'     => \$vcf_path, # path to vcftools vcf-compare script
                         'studyname=s'       => \$studyname);

if (!$list || !$studyname) {die $USAGE;}
$datadir ||=".";
$vcf_compare = $vcf_path if $vcf_path;

#===========================================================
# If tabix is not in the PATH, add its location to $PATH
#===========================================================
my $PATH = `echo \$PATH`;
my $reset_path;
my $tabix_check = `which tabix`;

if (!$tabix_check) {
 print STDERR "tabix not found, adding [$path_to_tabix] to the PATH...\n" if DEBUG;
 $ENV{PATH} = "$path_to_tabix:$PATH";
 $reset_path = 1;
} else {
 print STDERR "Found tabix, proceeding...\n" if DEBUG;
}


# ==========================================================
# First, if we have exisiting matrix, read values from there
# ==========================================================

my @old_matrices;

if (defined $oldmatrix && $oldmatrix=~/\S+/) {
 @old_matrices = grep{/\S+/} split(",",$oldmatrix);
}
 
if (@old_matrices && @old_matrices > 0) {
print STDERR "Found ".scalar(@old_matrices)." old matrices [$oldmatrix]\n" if DEBUG;
my $interact = 0;

foreach my $old (@old_matrices) {
 my %id_added = ();

 if ($old && -e $old) {
  print STDERR "We have data in ".$old." , will append our results to the existing matrix\n" if DEBUG;
  open(OLD,"<$old") or die "Couldn't read from the file with previously calculated indexes";
  my $sindex; # index of number of snps
  my $idstring = <OLD>;
  $idstring=~s/^0*//; # Funny entries should be removed by this command
  chomp($idstring);
  my @tempids = ();
  @tempids = grep {/\S+/} split("\t",$idstring);

  if ($tempids[$#tempids] =~ /SNP/) {
   $sindex = $#tempids;
   pop @tempids;
  }
  print STDERR "Header has ".scalar(@tempids)." ids\n" if DEBUG;
  map {&id_file($_)} (@tempids);
  map {$id_added{$_}++} (@tempids);

  while (<OLD>) {
   chomp;
   my @temp = split("\t");
   &id_file($temp[0]);
   $id_added{$temp[0]}++;
   $temp[0] = $` if $temp[0] =~/.snps.raw.vcf.gz$/;
   $temp[0] =~s!.*/!!;
   next if !$files{$temp[0]};
   my $idx = 1;
   T:
   foreach my $t (@tempids) {
     $t =~s!.*/!!;
     $t =~s!\..*!!;
     if (!$files{$t}) {next T;} # We should get rid of all funny entries here (like random 00 appering in files randomly)
     $t = $` if $t=~/.snps.raw.vcf.gz$/;
     $interact++ if ($temp[0] ne $t);
     $seen{$files{$temp[0]}}->{$files{$t}}++;
     $seen{$files{$t}}->{$files{$temp[0]}}++;

     $matrix{$files{$t}}->{$files{$temp[0]}} = $matrix{$files{$temp[0]}}->{$files{$t}} = $temp[$idx++] if ($files{$t} && $files{$temp[0]});
    }
   
   $snps{$files{$temp[0]}} = $sindex ? $temp[$#temp] : &calculate_snps($temp[0]);
  }
  print STDERR "Identified ".$interact." interactions using old matrices\n" if DEBUG;
  close OLD;
  print STDERR "We saw ".scalar(keys %id_added)." ids in this matrix\n" if DEBUG;
  print STDERR "We have ".scalar(keys %matrix)." unique ids in similarity matrix now...\n" if DEBUG;
  print STDERR "We have ".scalar(keys %files)." Files identified\n" if DEBUG;
 } else {
   next;
 }
}
print STDERR "Identified ".$interact." total interactions using old matrices\n" if DEBUG;
}

# =================================================================================
# Collect information on files, build matrix using new (and old, if available) data
# =================================================================================

if ($list=~/\,/) {
 # We have a comma-delimited list of files
 my @files = split(",",$list);
 map{&id_file($_)} @files; 
 $sublists[0] = [];
 map {push (@{$sublists[0]}, $_)} (keys %ids);
 $sublists[1] = $sublists[0];
} else {
 open(LIST,"<$list") or die "Cannot read from the list file [$list]";
 # make two sublists, depending on presence of empty line in the middle
 my $subidx = 0;
 $sublists[0] = [];
 $sublists[1] = [];
 
 while(<LIST>) {
  chomp;
 
  # Split the list if there's an empty line
  if (/^$/) {
   $subidx = 1;
   next;
  }
  push(@{$sublists[$subidx]}, &id_file($_));
 }

 if (!defined $sublists[1] || scalar(@{$sublists[1]}) == 0) {
   print STDERR "List is not split, will synchronize sublists\n" if DEBUG;
   $sublists[1] = $sublists[0];
 } 

 print STDERR "Got ".scalar(@{$sublists[0]})." and ".scalar(@{$sublists[1]})." elements in sublists\n" if DEBUG;
 close LIST;
}

print STDERR scalar(keys %ids)." genotypes collected\n" if DEBUG;

my $count = 1;

foreach my $id(@{$sublists[0]}) { #keys %ids) {
 # Take care of SNP number calculation
  
 my $sample = $ids{$id};

 print STDERR "Working on ".$count++." of ".scalar(@{$sublists[1]})." samples\n";
 SM:
 foreach my $s(@{$sublists[1]}) { 
   if (!$seen{$s} && !$oldmatrix) { # If we saw it, don't calculate snps
    $snps{$s} ||= &calculate_snps($ids{$s});
   }
   if ($s eq $id) {
    $matrix{$id}->{$s} = 1;
    next SM;
  }
  
  my $file1 = $datadir.$ids{$id};
  my $file2 = $datadir.$ids{$s};

  if ($seen{$id}->{$s}) {next;} 
  print STDERR "Interaction of $id and $s unseen, calculating...\n" if DEBUG;
  $seen{$id}->{$s}++;
  $seen{$s}->{$id}++;
  
  $file1 .= ".snps.raw.vcf.gz" if $file1 !~/gz$/;
  $file2 .= ".snps.raw.vcf.gz" if $file2 !~/gz$/;
  print STDERR "Will check files $file1 and $file2\n" if DEBUG;
  if (! -e $file1 || ! -e $file2){print STDERR "File(s) $file1 or $file2 not FOUND!\n";next;}
  print STDERR "Will run $vcf_compare $file1 $file2 ...\n" if DEBUG;

  # Touch .tbi files to make sure that index is newer than vcf
  `touch -h $file1.tbi || echo "Couldn't touch $file1.tbi"`;
  `touch -h $file2.tbi || echo "Couldn't touch $file2.tbi"`;

  my @compares = `$vcf_compare $file1 $file2 | grep \"^VN\"`;
  my @numbers = (0,0,0);

  for (my $i = 0; $i < @compares; $i++) {
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

my @heads = ();
map {/.snps.raw.vcf.gz/ ? push(@heads,$`) : push(@heads,$ids{$_})} (sort @{$sublists[0]}); 
print join("\t",("",@heads,"SNPs")); 
print "\n";

 foreach my $sample(sort @{$sublists[1]}) { 
 print $sample=~/.snps.raw.vcf.gz/ ? $` : $ids{$sample};
 TF:
 foreach my $ss(sort @{$sublists[0]}) { 
   my $value = $matrix{$sample}->{$ss} || $matrix{$ss}->{$sample};
   print $value ? "\t$value" : "\t0";
 }
 # Colors will be assigned using another script
 print $snps{$sample} ? "\t$snps{$sample}" : "\t0";
 print "\n";
 }


$ENV{PATH} = $PATH if $reset_path;

# ==================================================
# Subroutine for processing (registering) a vcf file
# ==================================================

sub id_file {
 my $file = shift @_;
 return if $file=~/^0+$/;
 $file=~s!.*/!!; 
 $file=~s!\..*!!;
 $file = $` if $file =~/.snps.raw.vcf.gz$/;
 my $id;

 if ($file=~/(\d+)_($studyname.\d+)_/ || $file=~/(\d+)_([A-Z]+.\d+)_/) {
  $id = $2.$1;
  $ids{$id} = $file;
  $files{$file} = $id;
 } elsif ($file=~/^($studyname\_[0-1]\d+).(\S+)/ || $file=~/([A-Z]+.\d+)\.(\S+)/) {
  $id = join("_",($1,$2));
  $ids{$id} = $file;
  $files{$file} = $id;
 } else {
  $ids{$file} = $file;
  $files{$file} = $file;
 }
 return $id ? $id : $file;
}

# ==================================================
# Subroutine for calculating number of SNPs
# ==================================================

sub calculate_snps {
 my $file = shift @_;
 if ($file!~/$datadir/) {
  $file = $datadir.$file;
 }
 if ($file !~/.snps.raw.vcf.gz$/) {
   $file.=".snps.raw.vcf.gz";
 }
 print STDERR "Calculating SNPs for $file\n" if DEBUG;

 my $result = `zcat $file | grep -v ^# | wc -l`;
 chomp($result);
 print STDERR $result." SNPs found\n" if DEBUG;
 return $result;
}
