#!/usr/bin/perl -w

# ==========================================================
# Script for sending alerts to a list of watcher
# Looks for 'FLAGGED' text in the alt tag for heatmap images
# in the report .html
# ==========================================================

use strict;
use Getopt::Long;
use constant DEBUG=>0;

my($html_path,$watchers,$bundle,$studyname); # external parameters
my(@affected_samples,@emails,$vetted_emails); # internal variables

my $MAIL = 'mail';
my $FLAGGED = '\"color:\#f93b08;\"';
my $MESSAGE = 
" Hi, it appears that  there is a possible sample swap detected for study STUDY,
Affected samples are:\nAFFECTED\nAnd the path to the report bundle is BUNDLE
please take action,\n\nSample Fingerprinting Workflow";

my $USAGE="plotReporter.pl --report=[req] --watchers=[req] --bundle=[req] --studyname=[req]";

my $result = GetOptions ('report=s'    => \$html_path,  # path to the report .html with HTML code to check
                         'watchers=s'  => \$watchers,   # comma-delimited list of emails of watchers
                         'bundle=s'    => \$bundle,     # final path to the report bundle
                         'studyname=s' => \$studyname); # We should include this into the body of email message

if (!$html_path || !$watchers || !$bundle || !$studyname) {die "Couldn't find valid set of arguments, $USAGE";}
chomp($studyname);
# =====================================================================
# 1 Validate emails, if we don't have a good list of vetted emails, die
# =====================================================================
@emails = split(",",$watchers);
$vetted_emails = join(",",&vet_emails(\@emails));
if (!$vetted_emails && $vetted_emails !~/\@/) {die "We don't have any valid email(s) for sending an alert to";}


# 2 Open report and detect FLAGGED Heatmaps=======================================================
# We rely on the assumption that by this time we nave an unvetted list of emails that is not empty
# At the end of this step , we should have an array of titles that correspond to flagged hetmaps
# ================================================================================================
open (HTML,"<$html_path") or die "Couldn't read from report [$html_path]";
@affected_samples = ();
while (<HTML>) {
 chomp;
 while (/label style=$FLAGGED>(\S+?)\<\/label/g) { 
    push (@affected_samples,$1);
 }
}
close HTML;

if (scalar(@affected_samples == 0)) {print STDERR "No flagged data found, stopping...\n";
                                     exit;}
my $affected = join(",\n",@affected_samples);
print STDERR "Affected Samples: \n".$affected."\n" if DEBUG;
$MESSAGE=~s/STUDY/$studyname/;
$MESSAGE=~s/AFFECTED/$affected/;
$MESSAGE=~s/BUNDLE/$bundle/;

print STDERR $MESSAGE;

# ====================================================
# We have a vetted list of recepients, send the email!
# ====================================================

`echo \"$MESSAGE\" | $MAIL -s "Possible Sample Swap/Mix-up Alert, Detected by Sample FIngerprinting" $vetted_emails`;

# Vetting emails
sub vet_emails {
 my $emails = shift @_;
 my @vetted = ();
 if (!$emails || @{$emails} == 0) {
  return undef;
 }
 return map {$_ if (/\w+\@.+/)} @emails;
}
