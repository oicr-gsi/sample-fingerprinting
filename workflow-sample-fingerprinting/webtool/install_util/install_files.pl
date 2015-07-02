#!/usr/bin/perl -w

use strict;
use Cwd;
use FindBin '$Bin';

my $origdir = cwd;
my $homedir = "$Bin/..";

chdir $homedir or die "couldn't cd to $homedir: $!\n";
my $basedir = shift @ARGV;
my @dirs    = ("cgi-bin","html");

print STDERR "Starting installing into $basedir\n";
chop($basedir) if $basedir =~m!/$!;

# Make directory tree
`mkdir -p $basedir/html/sampleswap/images`;
`mkdir -p $basedir/html/sampleswap/js`;
`mkdir -p $basedir/html/sampleswap/reports`;
`mkdir -p $basedir/cgi-bin/sampleswap`;

# Copy all files
my $wd = `pwd`;
chomp($wd);
print STDERR "Working dir is $wd\n";

&copy_dir("$wd/html",$basedir.'/html/sampleswap/',".html"); 
&copy_dir("$wd/html/images",$basedir.'/html/sampleswap/images/',".png");
&copy_dir("$wd/cgi-bin",$basedir.'/cgi-bin/sampleswap/');
&copy_dir("$wd/lib",$basedir.'/cgi-bin/sampleswap/',".pm");

system("cp","$wd/html/js/tool_scripts.js",$basedir.'/html/sampleswap/js/');

print STDERR <<END

***********************************************************************
Don't forget to do sudo chown www-data $basedir/html/sampleswap/reports
***********************************************************************
END
;

# ===========================
# copy a dir content
# ===========================
sub copy_dir {
 my($src,$dest,$wildcard) = @_;
 if (-d $src && -d $dest) {
   opendir(DIR,$src) or die "Couldn't read from directory [$src]";
   my @files = $wildcard ? grep{/$wildcard$/} readdir DIR : grep{!/\.$|\.\.$/} readdir DIR;
   close DIR;

   foreach my $file (@files) {
     `cp $src/$file $dest/`;
   }
 } else {
   warn "Something wrong with the directories, cannot copy from $src";
 }
}
