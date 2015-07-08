#!/usr/bin/perl 

#
# This is a development version of webtool that creates a custom report for Sample Fingerprinting workflow
#

use strict;
use warnings;
use CGI qw(:standard escape form start_table end_table);
use CGI::Carp qw(fatalsToBrowser);    # Remove for production use
use SwapReporter;

$CGI::POST_MAX = 1024 * 15000;  # maximum upload filesize is 15M
$CGI::DISABLE_UPLOADS = 0; # 1 disables uploads, 0 enables uploads
use constant SAMPLESPERSLICE=>8;

# Branch here if have either file and datadir or path to tempdir or nothing

if (my $file = param('upmatrix')) {
 my $dir = param('datadir') ;
 &handle_fileUpload($file,$dir);
} elsif (my $donors = param('donors')) {
 my $matrix  = param('matrix');
 my $repname = param('out_file');
 &generate_report($matrix,$donors,$repname);
} else {
 &draw_uploadForm();
}

exit 0;

# ========================================================
# File upload form 
# ========================================================
sub draw_uploadForm {
 my $datadir = "report_".$$;
 print header;
 print start_html;
 print h2("Download similarity matrix to start using Custom Report Tool"),
       br,
    table({-border=>"0", -cellspacing=>"10", -cellpadding=>"10", -bgcolor=>"lightgray"},
             TR({-align=>"LEFT", -valign=>"BOTTOM"},
                 td({-valign=>"top"},
                    form({-enctype=>"multipart/form-data", -method=>"POST"},
                         "Similarity Matrix to upload:",
                         input({-type=>"file",-name=>"upmatrix"}),br,br,
                         input({-type=>"hidden",-name=>"datadir",-value=>$datadir}),
                         submit({-name=>"Upload"}),
                         "(Press to upload the file)"
                         )
                   )
               )
       );
 print end_html;
}
# ========================================================
# Handle Upload
# ========================================================
sub handle_fileUpload {
 my ($file,$datadir) = @_;

 my $newdir = `pwd`;
 $newdir=~s/cgi-bin/html/;
 $newdir=~s/\s*$//;
 $newdir.="/reports/".$datadir;
 `mkdir -p $newdir`;
 $newdir.="/" if $newdir!~m!/$!;
 
 my $upfile = $file;
 my $basename = GetBasename($upfile);
 my $localfile = join("/",($newdir,$basename)); 
 my %donors = ();
 

 # give some feedback to browser
 warn "Saving the file $basename to $localfile\n";
 $datadir = $1 if $localfile=~m!^(\S+/)!;
 $datadir ||="./";
 
 # Handle compressed files here
 if ($upfile=~/\.zip$|\.tgz$|\.tar\.gz$/) {
  warn "Got a compressed file, going to extract data\n";
  if (! open(ZIPFILE, ">$localfile") ) {
	warn "Can't open $localfile for writing - $!";
	exit(-1);
 }
 # Transfer it first  
  my $nBytes = 0;
  my $totBytes = 0;
  my $buffer = "";

  binmode($upfile);

  while ( $nBytes = read($upfile, $buffer, 1024) ) {
 #while ( $nBytes = read($fh, $buffer, 1024) ) {
	print ZIPFILE $buffer;
	$totBytes += $nBytes;
 }

 close(ZIPFILE);
  
  if ($localfile=~/\.zip$/) {
   `unzip $localfile -d $datadir`;
  } elsif ($localfile=~/\.tgz$|\.tar\.gz$/) {
  `tar -C $datadir -xzf $localfile`;
  } else {
   die "Could not unzip the input bundle!";
  }

  # find the matrix file and set self->matrix to that value
  warn "Looking for similarity matrix in $datadir";
  my @files = grep {/csv|txt/} `find $datadir -name \*jaccard.matrix\*`; #`ls --color="none" $datadir\*matrix\* | perl -ne '{chomp;s/.* //g;print \$_}'`;
  map{chomp($_)} @files; 
  if ($files[0] =~ /matrix/) {
   $localfile = $files[0];
   warn "Got matrix file $localfile";
  } else {
   die "Couldn't find any similarity matrix files in your input";
  }  
  
 } else {  
   if (! open(OUTFILE, ">$localfile") ) {
	warn "Can't open $localfile for writing - $!";
	exit(-1);
   }
   
   while (<$upfile>) {
 	s/\r//;
	print OUTFILE "$_";
	}
	close OUTFILE;
}

 
 my $first = 1;
 my $fields;
 if (! open(OUTFILE, "<$localfile") ) {
	warn "Can't open $localfile for reading - $!";
	exit(-1);
 }
 
 while (<OUTFILE>) {
 	# Do not extract donors from the first line	
	if ($first) {
	 $first = 0;
	 next;
	}

   # Check the number of fields
   my @temp = split("\t");
   $fields ||= scalar(@temp);
   
   if ($fields != scalar(@temp)) {
     print header;
     print start_html,
           h2("There was an error loading file $basename, number of field is not the same for all lines"),
           end_html;
     exit(-1);
   }

   # All is well, extract the donor name
	if (/^SWID_\d+_([A-Z]+_\d+)_/ || /^([A-Z]+_\d+)_/ ) {
     $donors{$1} = 1;	
	}
 }
 
 close OUTFILE;
 my $list = join(",",(keys %donors)); 
 $file = $localfile;
 $file =~s!$datadir!!;

 &show_customTool($list,"options",$newdir,$file);
 
}


# ========================================================
# Custom tool (with redirect to report page)
# ========================================================
sub show_customTool {

 my($list,$options,$dir,$matrix) = @_;
 my @donors = grep{/\S+/} split(",",$list);
 return if @donors == 0;
 
 print header;
 print start_html(-title=>'Custom Report Creation',
                  -author=>'pruzanov@oicr.on.ca',
                  -meta=>{'keywords'=>'sample swap detection genotype fingerprinting',
                                      'copyright'=>'&copy; 2013 OICR'},
                  -script=>[{-type => 'text/javascript',
                             -src  => '../../html/sampleswap/js/tool_scripts.js'}],
                  -BGCOLOR=>'white');

 print h2("Choose Donors for custom report");
 print br;
 my @inputs = ();

 for (my $d=0; $d<@donors; $d++) {
  push(@inputs,(checkbox({-id=>'check'.$d,
                         -name=>"$donors[$d]",
                         -checked=>0,
                         -value=>"$donors[$d]",
                         -onclick=>"countDonors(".scalar(@donors).",\'$options\',".SAMPLESPERSLICE.")",
                         -height=>"32"}),br));
 }
 
 my @tool_tds = (td({-valign=>"top",-width=>"150"},@inputs),
                 td({-valign=>"top"},
                   "Note that only maximum of 8 donors can be chosen to generate report",                   
                   br,
                 form({-name=>"genrep",-enctype=>"multipart/form-data", -method=>"POST"},
                   hidden({-name=>"donors",-value=>"Donors"}),
                   hidden({-name=>"matrix",-value=>join("/",($dir,$matrix))}),
                   hidden({-name=>"out_file",-value=>"custom_report_1234"}),                 
                   br,br,
                   submit({-name=>"Generate"}),
                   br,br,
                   "(After pressing Generate an HTML report will be created and shown here)"
                  )));
                   
 print table({-border=>0,-cellspacing=>10,-cellpadding=>10,-bgcolor=>"lightgray"},
              Tr({-align=>'LEFT',-valign=>'BOTTOM'},@tool_tds));

 print end_html;
}

# ===================================================================
# Generate report using parameters passed by customization interface
# ===================================================================

sub generate_report {
 my($matrix,$donors,$repname) = @_;
 my $sw = new SwapReporter($matrix,$donors,$repname);
 $sw->generateReport();
}


# ========================================================
# GetBasename - delivers filename portion of a fullpath.
# ========================================================

sub GetBasename {
	my $basename = shift;
   
	# check which way our slashes go.
	if ( $basename =~ /(\\)/ ) {
		$basename=~s!.*\\!!;
	} else {
		$basename=~s!.*/!!;
	}
	return $basename;
}


