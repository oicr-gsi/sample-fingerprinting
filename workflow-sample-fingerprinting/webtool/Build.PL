#!/usr/bin/perl -w

use strict;
use FindBin '$Bin';

use lib "$Bin/install_util";
use File::Copy;
use SwapReportInstall;

  my $build = SwapReportInstall->new
    (module_name => 'SwapReporter',
     license  => 'perl',
     dist_author        => 'Peter Ruzanov <pruzanov@gmail.com>',
     dist_abstract      => 'SwapReporter webtool for SampleFingerprinting workflow at OICR',
     requires => {
                  'perl'          => '5.10.0',
                  'CGI'  => '1.23'
                 }
    );
  $build->create_build_script;


print STDERR <<END;

**** WARNING *****************************************************************************************
This is rather unflexible setup, it requires user to specify a top-level directory for installation,
usually \$HOME/public_html or /var/www. The module and all the scripts with .html documents and images
will be installed in directory specified at the installation time. Changing ownership for the reports dir
will be required (apache needs to write into this directory)
ALSO, R needs to be installed! R scripts won't run if there's no Rscript in the PATH
******************************************************************************************************
END
