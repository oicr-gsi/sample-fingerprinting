#!/usr/bin/perl -w

package SwapReporter;

our $VERSION="1.0a";

use strict;
use warnings;
use CGI qw/:standard form/;
use IO::File;
use Data::Dumper;
use constant DEBUG=>0;
use constant WEBRUN=>1; # if used on command-line this should be set to 0
use constant THRESHOLD=>30; # That many SNPs every genotype should have
use constant SAMPLESPERSLICE=>8;
use constant PNGSIZE=>750;

=head2 SUMMARY
 
 SwapReporter module works with similarity matrix and a list of donors. Produces html report with heatmap images

=head2 USAGE
 

 my $swr = new SwapReporter(matrix_file,$donors);
 $swr->generateReport();

 Replicates make_report script and may automatically split data into smaller chunks, produce multiple heatmaps etc.
 Optionally it uses tempdir parameter which may point to a directory with coverage depth data if that is different than datadir
 (web version should not create fingerprints if there's no finfiles dir)

=head2 EXAMPLES

 my $swr = new SwapReporter(matrix_file,$donors);
 $swr->generateReport();

 my $swr = new SwapReporter("PCSI_WG.matrix.all.txt","PCSI_0001,PCSI_0002,PCSI_0003");

 $swr->generateReport();

 generateReport() will produce all heatmap images, html report, csv files with genotype info and similarity info, fingerprint images (if this is not web version)

=cut


sub new {
    my $class = shift;
    my @args = @_;
    my @dons = (split ",",$args[1]);
    my %donors = map{$_=>1} @dons;
    my $self = bless { matrix    => $args[0],
                       fmatrix   => "",         # filtered matrix file
                       donors    => \%donors,
                     # repname   => $args[2],
                       reports   => {},
                       samples   => {},
                       ids       => {},
                       snpindex  => 0,
                       studyname => "",
                       datadir   => "./",
                       filtered  => {},
                       finfiles  => 0}, ref $class || $class;    # will set finfiles flag to 1 if we have it (and if we do, genotype file and fingerprint popup will be created)
}

 


sub generateReport {
 my $self = shift;
 $self->loadInput();
 $self->loadMatrix();
 $self->createHeatmap();
 $self->createHtml();
}

# ================================================================================================================
# Setters/Gettersi TODO: make one for ids, samples etc
# ================================================================================================================
sub finfiles {
 my($self,$flag) = @_;
 $flag ? $self->{finfiles} = $flag : return $self->{finfiles};
}

# ================================================================================================================
# Loading input either as zipped bundle (matrix file and fin files in finfiles directory) or just matrix
# ================================================================================================================
sub loadInput {
 my $self = shift;

 # Find out what we have
 print STDERR "Matrix is $self->{matrix}\n" if DEBUG;
 my $dir = $1 if $self->{matrix}=~m!^(\S+/)!;
 $dir ||="./";
 $self->{datadir} = $dir;
 
 # find directories finfiles and fingerprints and check that they have .fin and .png files
 my $numfin = `find $dir -name \*\.fin | wc -l`;
 chomp($numfin);
 if ($numfin > 0) {$self->finfiles(1);}

}

# =================================================================================================================
# Script for filtering and clustering genotypes, splitting matrix into a set of smaller chunks and producing images 
# + html wrapper
# =================================================================================================================
# Modularized load function
sub loadMatrix {
 my $self = shift @_;
 my @lines;

 if ($self->{matrix} && -e $self->{matrix}) {
 my $fh = new IO::File("<$self->{matrix}") or $self->printError("There was an error reading file [".$self->{matrix}."]");
 
 my $firstline = <$fh>;
 chomp($firstline);
 my @heads = split("\t",$firstline);
 map {if ($heads[$_]=~/^SNP/){$self->{snpindex} = $_}} (0..$#heads);

 $self->{snpindex} or $self->printError("The matrix file [".$self->{matrix}."] is missing the column with number of SNPs called per genotype");
 @lines = ($firstline);

 # Find what is the index of 'SNPs' column and use it later to remove those entries with fewer than THRESHOLD SNPs
 while (<$fh>) {
  #chomp;
  my @temp = split("\t");
  my $trimmed_name = $temp[0];
  $trimmed_name=~s/(SWID_\d+)_(.*)/$2\_$1/;
  my($don,$studyname);
  if ($trimmed_name=~/^((\S+?)_\d+)/) {
    $self->{studyname} ||= $2;
    $don = $1;
   
  } 

  # TODO If there are no donors, we should not use self->{donors} in filtering

  if ($temp[$self->{snpindex}] && $temp[$self->{snpindex}] >= THRESHOLD && $don && $self->{donors}->{$don} && $temp[0]=~/(\d+)_$don/) {
    $self->{ids}->{$temp[0]} = join("_",($don,$1));
    $self->{samples}->{$self->{ids}->{$temp[0]}} = {sample=>$don,
                                                    file  =>$temp[0],
                                                    name  =>$trimmed_name}; # register a file as pertaining to a certain sample (studyname_sampleid)
    push(@lines,join("\t",@temp));
  } else {
    $self->{filtered}->{$temp[0]} = $trimmed_name;
    next;
  }
 } # reading from matrix ends here

 $fh->close;

 # Post-filtering: register filtered files as skipped
 # Open a file for writing out filtered matrix (with friendlier file ids):
 my $mfile = $self->{datadir}."matrix_filtered_$$.csv";
 $fh->open(">$mfile") or $self->printError("Couldn't write filtered data into a file [$mfile]");

 my @filterhead;
 HEAD:
 foreach my $hf (@heads) {
  if ($self->{filtered}->{$hf}) {next;}
  push(@filterhead,$hf) if $hf=~/\w+/;
 }

 print STDERR "Filtered header retained ".scalar(@filterhead)." of ".scalar(@heads)." original ones, skipped indexes: ".scalar(keys %{$self->{filtered}})."\n" if DEBUG;

 map{print $fh "\t".$_ if $self->{ids}->{$_}} (@filterhead);
 print $fh "\tSNPs\n";

 # Get read of filtered files (lines in the matrix)
 LINE:
 foreach my $line (@lines) {
  my @tlines=split("\t",$line);
  if ($line=~/^(\S+)\t/) {
    if ($self->{filtered}->{$tlines[0]}){next LINE;}

    print $fh $tlines[0];
    IDX:
    foreach my $line_idx (1..$#tlines) {
      if ($self->{filtered}->{$heads[$line_idx]}) {next IDX;}
      $tlines[$line_idx] =~/NA/ ? print $fh "\t0" : print $fh "\t$tlines[$line_idx]";
    }
  } else {
    next LINE;
  }
 }
  $fh->close;
  $self->{fmatrix} = $mfile;
 } else {
  $self->printError("No valid matrix file supplied, I cannot continue with no input");
 } 
}

# ==========================================================
# Load file (filtered matrix), format data and make heatmaps
# ==========================================================
sub createHeatmap {
 my $self = shift;
 my $matrix = $self->{fmatrix};
 my %colors;
 my @colors = qw/red orange yellow green lightblue blue purple darkgreen black pink/;

 # Just use the file here, make heatmap 
 my $pngfile = $';
 $pngfile ||=$matrix;

 my $matfile  = $pngfile;
 $pngfile =~s/csv$/png/;
 $matfile =~s/filtered/formatted/;

 # Temporary matrix file for a slice
 my $fm = new IO::File("<$matrix") or die "Cannot read from file [$matrix]";
 my $fo = new IO::File(">$matfile") or die "Cannot write to file [$matfile]";
 my $first = <$fm>;
 my @names = split "\t",$first;
 shift @names; # remove 1st (useless) element

 my %indexes;
 my $colcount = 0;

 NAME:
 for (my $i = 0; $i < @names; $i++) {
 if ($self->{ids}->{$names[$i]}) {
  $names[$i]=~s!.*/!!;
  print $fo "\t".$self->{samples}->{$self->{ids}->{$names[$i]}}->{name};
  $indexes{$i} = $self->{ids}->{$names[$i]};
  $colors{$self->{samples}->{$self->{ids}->{$names[$i]}}->{sample}} ||= $colors[$colcount++];
  }
 }

 print STDERR "Assigned ".scalar(keys %colors)." colors to clusters\n" if DEBUG;
 print $fo "\tSNPs\tColor\n";

 while(<$fm>) {
  chomp;
  my @temp = split("\t");
  $temp[0]=~s!.*/!!; # remove path, leave the name
  my $col = $colors{$self->{samples}->{$self->{ids}->{$temp[0]}}->{sample}};
  print $fo $self->{samples}->{$self->{ids}->{$temp[0]}}->{name};
  
  foreach my $idx(sort {$a<=>$b} keys %indexes) {
   if ($temp[$idx + 1] ne "NA") {
      print $fo "\t".$temp[$idx + 1];
   } else {
      print $fo "\t0";
   }
  }

  print $fo "\t".join("\t",($temp[$#temp],$col)),"\n";
 }
 $fo->close;
 $fm->close;

 # Produce images
 my $flagged = 'FALSE';
 my $title = join(",",(keys %{$self->{donors}}));

 print STDERR "Will Rscript create_heatmap.r $matfile \"$title\" 400 $pngfile ".PNGSIZE." $flagged\n" if DEBUG;
 my $size = PNGSIZE;
 my $clustered_ids =  `Rscript create_heatmap.r $matfile \"$title\" 400 $pngfile $size $flagged`;
 my @clustered_ids = grep {/$self->{studyname}/} split(" ",$clustered_ids); 

 my @fingers = ();

 for (my $cl = 0; $cl < @clustered_ids; $cl++) { 
  ID:
  foreach my $id (keys %{$self->{ids}}) {
    if ($clustered_ids[$cl] eq $self->{samples}->{$self->{ids}->{$id}}->{name}) {

    # Depending on wheather we have finfiles dir or not, create fingerprint images

     if ($self->finfiles()) {
     print STDERR "Creating $cl of ".scalar(@clustered_ids)." fingerprint images\n" if DEBUG;
     my $png = $self->{datadir}.$self->{ids}->{$id}.".fp.".$cl.".png";
     
     my $fin = $id.".fin";
     $fin =~s!.*/!!;

     print STDERR "Will Rscript create_fingerprints.r $self->{datadir}finfiles/ $fin $colors{$self->{samples}->{$self->{ids}->{$id}}->{sample}} 400 $png \n" if DEBUG;
     `Rscript create_fingerprints.r $self->{datadir}finfiles/ $fin $colors{$self->{samples}->{$self->{ids}->{$id}}->{sample}} 400 $png`;
     $png=~s!.*(report_\d+)/!!;
     $png=~s!.*/!!;
     push(@fingers,{img=>$png,
                    name=>$self->{samples}->{$self->{ids}->{$id}}->{name}});
     last ID;
     }
    }
  }
 }

 # Register the image name in the report hash
 $pngfile=~s!.*(report_\d+)!../../html/sampleswap/reports/$1!;
 $matfile=~s!.*(report_\d+)!../../html/sampleswap/reports/$1!;
 $self->{reports} = {img=>$pngfile,
                     fp=>[@fingers],
                     flagged=>$flagged eq "TRUE" ? "FLAGGED" : "OK",
                     matrix=>$matfile,
                     title=>$title};
 if ($self->finfiles()) {
    $self->printout_snps();
 }
 print STDERR Dumper($self->{reports}) if DEBUG;
}

# ============================================================================
# Multi-heatmap code gets removed from here temporary, to be implemented later
# Using R (heatmap) cluster samples, 
# we need the re-arranged list for next step
# ============================================================================

# =================================================================
# make HTML report (will call a couple of subroutines)
# =================================================================

# These images are hardcoded, not supposed to be customizable

sub createHtml {
 my $self = shift;
 my $helplink = "../../html/sampleswap/help.html";
 print header;
 print start_html(-title=>'Sample Fingerprinting Report',
                  -author=>'pruzanov@oicr.on.ca',
                  -meta=>{'keywords'=>'sample swap detection genotype fingerprinting',
                                      'copyright'=>'&copy; 2013 OICR'},
                  -script=>[{-type => 'text/javascript',
                            -code  => 'function showFingerprints(snapshot){window.open(snapshot,"_blank","width='.(PNGSIZE + 1).',height=600,toolbar=0,menubar=0,status=1,scrollbars=yes,resizable=1")}'}],
                  -BGCOLOR=>'white');
 print button(-onClick=>"showFingerprints('$helplink')",
              -name=>"help_button",
              -value=>"Help");
 print "\n&nbsp;&nbsp;\n";
 my $matrix_link = WEBRUN ? $self->{fmatrix} : $self->{matrix};
 $matrix_link=~s!.*(report_\d+)!../../html/sampleswap/reports/$1!;
 print button(-onClick=>"window.location.href=\'$matrix_link\'",
              -name=>"download_button",
              -value=>"Download Data");
 print h2("Sample Fingerprinting for ".$self->{studyname}." study");
 print br;
 #1. Suspicious samples
 if (scalar(keys %{$self->{flagged}->{files}}) > 0) {
  my @flagged = map{Tr({-align=>'LEFT',-valign=>'BOTTOM'},td($_))} (keys %{$self->{flagged}->{files}});

  print h3("Files flagged as potential sample swaps:");
  print br;
  print table({-border=>0},
               @flagged);
}


 #2. Filtered files:
 if (scalar(keys %{$self->{filtered}}) > 0) {
  my @filtered = map{Tr({-align=>'LEFT',-valign=>'BOTTOM'},td($_))} (values %{$self->{filtered}});

  print h3("Files skipped due to low coverage/small number of SNPs or single file in a sample:");
   print br;
   print table({-border=>0},
                @filtered);
 }


 # Define table = 3 columns always, rows - depending on the number of heatmaps
 my $n_rows = scalar(keys %{$self->{reports}})/3;
 $n_rows = int($n_rows) < $n_rows ? int($n_rows + 1) : int($n_rows);

 # 3. image of the heatmap and 4. button with a link to popup with fingerprints
 my @hmaps = $self->heatmap_rep($self->{reports});
 my @tab_rows = map {3*$_+2 <= $#hmaps ? Tr({-align=>'LEFT',-valign=>'TOP'},@hmaps[3*$_..3*$_+2])
                                      : Tr({-align=>'LEFT',-valign=>'TOP'},@hmaps[3*$_..$#hmaps]);} (0..$n_rows-1);

 print h3("Heatmaps based on similarity matricies:");
 print table({-border=>0},
            @tab_rows);

 print end_html;
}
#===========================================================================================
# process slices with two R scripts - one for heatmap, one for 'barcode'-looking fingerprint
# Here we are printing to two files, slightly different header and file ids
#===========================================================================================
# =================================================================
# Creates HTML for a table cell (heatmap, button for opening popup)
# =================================================================

sub heatmap_rep {
 my $self = shift;

 my $icondir = "../../html/sampleswap/images/";

 my $link_image = $icondir."fp_button.png";
 my $sim_image  = $icondir."sim_button.png";
 my $gen_image  = $icondir."gen_button.png";

 my $popup;
 if ($self->finfiles()) {$popup = $self->create_popup()};

 my @samples = split(",",$self->{reports}->{title});
 my @labels  = map{$self->{flagged}->{samples}->{$_} ? br.label({-style=>"color:#f93b08;"},$_) : br.label($_)} @samples;

 $self->finfiles() ? 
     return     td(img({-src=>$self->{reports}->{img},-width=>500,-height=>500,-alt=>"Heatmap_".$self->{reports}->{flagged}}),br,
          image_button({-src=>$link_image,-width=>111,-height=>32,-alt=>"clickglyph",-name=>"fingers",-onClick=>"showFingerprints(\'$popup\')",-value=>1,-align=>'MIDDLE'}),
                     a({-href=>$self->{reports}->{matrix},-target=>'_new'},img({-src=>$sim_image,-width=>111,-height=>32,-alt=>"matrix_glyph",-name=>"matrix",-value=>1,-align=>'LEFT'})),
                     a({-href=>$self->{reports}->{genotype},-target=>'_new'},img({-src=>$gen_image,-width=>111,-height=>32,-alt=>"genotype_glyph",-name=>"genotype",-value=>1,-align=>'LEFT'})),br,
                     @labels) :
     return     td(img({-src=>$self->{reports}->{img},-width=>500,-height=>500,-alt=>"Heatmap_".$self->{reports}->{flagged}}),br,           
                     a({-href=>$self->{reports}->{matrix},-target=>'_new'},img({-src=>$sim_image,-width=>111,-height=>32,-alt=>"matrix_glyph",-name=>"matrix",-value=>1,-align=>'LEFT'})),br,
                     @labels)
}

# ======================================================================================
# Function for writing HTML for a popup (for a cluster #, expects image being available)
# ======================================================================================

sub create_popup {
 my $self = shift;

 if (!$self->{reports}->{fp}) {
  return "../../html/sampleswap/error_popup.html";
 }

 my $popname = $self->{datadir}."fingerprints_popup_$$".".html";
 my $pop = new IO::File(">$popname") or die "Cannot write to [$popname]";

 print $pop start_html(-title=>"Individual Fingerprints Report",
                       -author=>'pruzanov@oicr.on.ca',
                       -meta=>{'keywords'=>'sample swap detection genotype fingerprinting',
                                           'copyright'=>'&copy; 2013 OICR'},
                       -BGCOLOR=>'white');

 my @t_rows = ();
 foreach my $finger (@{$self->{reports}->{fp}}) {
  my $nopath_name = $finger->{name};
  $nopath_name=~s!.*/!!;

  if ($self->{flagged}->{files}->{$finger->{name}}) {
   push(@t_rows,Tr({-align=>'LEFT', -valign=>'BOTTOM'},(td(img({-src=>$finger->{img},-alt=>'Open this report in the dir with the original report to see the fingerprint images'.$finger->{img}})),
                                                        td({'style'=>"font-size:small; color:#ff0000;"},label($nopath_name)))));
   print STDERR "FLAGGED File found!\n" if DEBUG;
  }else{
   push(@t_rows,Tr({-align=>'LEFT', -valign=>'BOTTOM'},(td(img({-src=>$finger->{img},-alt=>'Fingerprint '.$finger->{img}})),
                                                        td({'style'=>"font-size:small;"},label($nopath_name)))));
  } 
 }

 print $pop table({-border=>0,
                   -cellpadding=>0,
                   -cellspacing=>0},
                   @t_rows);

 print $pop end_html;
 $pop->close;
 $popname=~s!.*html/sampleswap!../../html/sampleswap!;
 return $popname;
}

#==========================================================================================================================
# A subroutine for printing out genotype report for a heatmap (slice) - will list all SNPs in checked 'hotspots' in a table
#==========================================================================================================================
sub printout_snps {
 my $self = shift;

 # Open .fin file and read info from there, build matrix for all files in the slice and print into a file
 my $fh_fin = new IO::File();
 my %snpinfo;
 my %snpcalls;
 my @titles = (); # Cell titles

 foreach my $id (sort keys %{$self->{ids}}) {
  my $finfile = $self->{datadir}."finfiles/".$id.".fin";
  my $file_ok = 1;
  $fh_fin->open($finfile) or $file_ok = 0;
  
  if (!$file_ok) {
   warn "File with snp info for $id is not available";
   next;
  }

  my $first = <$fh_fin>;
  if ($first!~/^CHROM/){next;}

  @titles = split("\t",$first) if !@titles;
  my $filename = $self->{samples}->{$self->{ids}->{$id}}->{file};
  while (<$fh_fin>) {
   chomp;
   my @temp = split("\t");
   $snpinfo{$temp[0]}->{$temp[1]} ||= $temp[2];
   $snpcalls{$temp[2]}->{$self->{samples}->{$self->{ids}->{$id}}->{file}} = $temp[3];
  }
  $fh_fin->close;
 }

 # Having collected all snp calls from the .fin files let's create a genotype report file
 my $fname = $self->{datadir}.join("_",($self->{studyname},"genotype_report_$$.csv"));
 $fh_fin->open(">$fname") or die "Couldn't write genotype report to [$fname]";
 print $fh_fin join("\t",@titles[0..2]);
 my @fnames = grep {/\S+/} map{$self->{samples}->{$_}->{name}} (sort keys %{$self->{samples}});
 print $fh_fin "\t",join("\t",@fnames),"\n";

 foreach my $chrom (sort keys %snpinfo) {
  foreach my $pos (sort {$a<=>$b} keys %{$snpinfo{$chrom}}) {
   print $fh_fin join("\t",($chrom,$pos,$snpinfo{$chrom}->{$pos}));
   FILE:
   foreach my $file_id (sort keys %{$self->{samples}}) {
    if ($self->{filtered}->{$self->{samples}->{$file_id}->{file}}) {next FILE;}
    $snpcalls{$snpinfo{$chrom}->{$pos}}->{$self->{samples}->{$file_id}->{file}} ? print $fh_fin "\t".$snpcalls{$snpinfo{$chrom}->{$pos}}->{$self->{samples}->{$file_id}->{file}} : print $fh_fin "\t";
   }
   print $fh_fin "\n";
  }
 }
 $fh_fin->close;
 $fname=~s!.*/(report_\d+?/)!../../html/sampleswap/reports/$1!;
 $self->{reports}->{genotype} = $fname;

}


# ========================================================================================
# If requested, we can provide the whole report in zipped archive
# ========================================================================================
sub createZip {
 print STDERR "Zipping is not implemented yet\n";
}

# ========================================================================================
# Print an error in HTML format
# ========================================================================================
sub printError {
 my $self = shift;
 my $message = shift;
 print header;
 print start_html;
 print h2($message);
 print end_html;
 exit 1;
}

1;
