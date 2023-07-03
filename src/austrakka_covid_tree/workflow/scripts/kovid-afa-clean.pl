#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std;
use File::Basename;
use File::Temp;
use Data::Dumper;

use FindBin;
use lib "$FindBin::RealBin/../perl5";

#......................................................................................
# globals


# used by kovid_enter/exit and tempdir
my $old_dir;

#....................................................................
sub msg { say STDERR "@_"; }

#....................................................................
sub err { msg("ERROR:", @_); exit(1); }

#....................................................................
sub wrn { msg("WARNING:", @_); }


sub cpus { 
  my($c) = qx(nproc --all);
  chomp $c;
  return $c || 1;
}

#....................................................................
sub fixfile { 
  my($f, $text) = @_;
  $text ||= "bummer";
  $f or err("$text - empty filename");
  -r $f or err("$text - can't read file '$f'");
  return abs_path($f);
}

#....................................................................
sub kovid_enter {
  msg("kovid toolkit by \@torstenseemann");
  msg(@_) if @_;
}

#....................................................................
sub kovid_exit {
  #chdir("/"); # allow tmpdir removal
  if ($old_dir) {
    msg("Changing back to dir: $old_dir");
    chdir $old_dir;
  }
  msg(@_) if @_;
  msg("Done.");
  exit(0);
}



my $EXE = basename($0);
my $IUPAC = "ACGTUWSMKRYBDHV";
my $DNA = "ACGT";

#......................................................................................
# defaults
my %opt = (
  'm' => 'N',
  'c' => 'upper',
);

#......................................................................................

sub usage {
  my($errcode) = @_;
  $errcode ||= 0;
  my $ofh = $errcode ? \*STDERR : \*STDOUT;
  print $ofh 
    "NAME\n  $EXE\n",
    "SYNOPSIS\n  Clean alignment chars\n",
    "USAGE\n  $EXE [options] all.afa > var.afa\n",
    "OPTIONS\n",
    "  -h       Print this help\n",
    "  -m CHAR  For missing data [$opt{m}]\n",
    "  -g       Keep gaps, don't covert to -m\n",
    "  -a       Keep ambig, don't convert to -m\n",
    "  -c CASE  Case: lower,upper,keep [$opt{c}]\n",
    "  -F       Put into 'FastTree' mode\n",
    "  -Q       Put into 'IQTree' mode\n",
    "END\n";
  exit($errcode);
}

#......................................................................................
# OPTIONS
getopts('hc:gam:FQ', \%opt) or exit(-1);
$opt{h} and usage(0);

if ($opt{'F'}) {
  msg("Using settings for FastTree");
  $opt{'m'} = '-';
  $opt{'a'} = undef;
  $opt{'g'} = undef;
}
elsif ($opt{'Q'}) {
  msg("Using settings for IQtree");
  $opt{'m'} = 'N';
  $opt{'a'} = 1;
  $opt{'g'} = 1;
}

$opt{c} or err("Need -c");
my $case = uc substr($opt{c}, 0, 1);
my $miss = $opt{'m'} or err("Need -m");
length($miss)==1 or err("-m $miss must be 1 char");

my $alpha = $opt{'a'} ? $IUPAC : $DNA;
my $keepgaps = $opt{'g'};
$alpha .= $opt{'m'};
# want '-' at end of regexp
$alpha .= '\-' if $keepgaps;

#......................................................................................

kovid_enter("This is $EXE");

msg("Case        :", $case);
msg("Missing is  :", $miss);
msg("Gaps become :", $keepgaps ? '-' : $miss);
msg("Alphabet    :", $alpha);

msg("Processing alignment...");
my $nseq=0;

while (<ARGV>) { 
  chomp;
  #print STDERR "#OLD $_\n";
  if (m/^>/) {
    $nseq++;
    #msg("DOne $nseq ...") if $nseq%1000==0;
  }
  else {
    # RNA to DNA
    s/U/T/gi;
    # gaps
    # this might be redundant now
    # as we put - in $alpha
    s/-/$miss/g unless $keepgaps;
    # alphabet
    s/[^$alpha]/$miss/gi;
    # case
    if ($case eq 'U') { $_ = uc($_) }
    elsif ($case eq 'L') { $_ = lc($_) }
  }
  print $_,"\n"; 
  #print STDERR "#NEW $_\n";
}
msg("Processed $nseq sequences.");

kovid_exit("Done.");

#.......................................
