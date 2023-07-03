#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use File::Which; 
use File::Temp;
use Cwd qw(abs_path);
use Getopt::Std;
use File::Copy;
use List::Util qw(sum min max);
use Cwd;
use FindBin;
use lib "$FindBin::RealBin/../perl5";

#......................................................................................
# globals
my $kovid_afa_clean = "$FindBin::Bin/kovid-afa-clean.pl";

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

#....................................................................
sub exe {
  my($tool) = @_;
  my $path = which($tool);
  $path ? msg("Found '$tool' - $path") 
        : err("Can't find '$tool' - please install it");
  return $path;
}


my $EXE = basename($0);

# this should be > log_10(29903)
my $PRECISION = 6;
my $BLMIN = 10**(-$PRECISION-1);
msg("BRLEN = $BLMIN");
# raxml-ng
#  --blmin VALUE 
# minimum branch length (default: 1e-6)
  
#......................................................................................
# defaults

my %opt = (
  'e'=>'FastTreeMP',
  'x'=>'',
  'F' => undef,  # use -fastest if enabled
  'B' => 1000,   # 0 = -nosupport 
  'p'=>'',
  'j'=> 8, #cpus(),
  'o'=> 'MN908947.3',
  'r'=> '../resources/MN908947.3.fna',
  't'=>'',
);

sub usage {
  print <<"EOHELP";
USAGE
  $EXE [options] -a in.afa -p out 
OPTIONS
  -a AFA    Input alignment file
  -p PREFIX Output file prefix
  -j CPUS   Threads to use [$opt{j}]
  -r FASTA  Reference FASTA [$opt{r}]
  -e EXE    FastTeee binary [$opt{e}]
  -x OPTS   Extra EXE opts [$opt{x}]
  -B NUM    Bootstraps [$opt{B}]
  -F        Faster but less accurate
  -t DIR    Working folder
EOHELP
  exit(0); 
}
getopts('ha:p:j:r:e:x:t:B:F', \%opt) or exit(-1);

$opt{h} and usage();
my $afa = fixfile($opt{a}, "need -a AFA");
my $ref = fixfile($opt{r}, "Need -r REF");
my $cpus = $opt{j} || 1;
my $out = $opt{p} or err("need -p PREFIX");
#my $outgroup = $opt{'o'};
my $ftree = $opt{'e'} or err("need -e EXE");

kovid_enter("This is $EXE");
exe($ftree);
exe("raxml-ng");
# exe($kovid_afa_clean);
exe($_) for qw(goalign gotree nw_stats nw_order nw_reroot);

my $refseq = fasta2array($ref);
my $L = length($refseq->[0][1]);
my $outgroup = $refseq->[0][0];
msg("Reference $outgroup is $L bp");

my $tdir = -d $opt{t} ? $opt{t} : tmpdir();
msg("Working folder: $tdir");
my $tmp = "$tdir/kovid.FT";
msg("Temp prefix:", $tmp);
#msg("CPUS:", $cpus);

msg("Purifying AFA for $opt{e}");
run("$kovid_afa_clean -F '$afa' > $tmp.afa1");

msg("De-duping AFA...");
run(
  "goalign dedup -i $tmp.afa1"
 ." -o $tmp.afa2"
);
msg("Unique seqs  :", num_lines("$tmp.dup")) ;

msg("Removing invariant sites...");
run(
  "goalign compress"
  ." -i $tmp.afa1 -o $tmp.afa3"
  ." --weight-out $tmp.afa3.weights"
);
my $var = num_lines("$tmp.afa3.weights");
msg("Unique sites : $var"), 

msg("Building tree with $ftree");

my $opts = $opt{x} // '';
$opts .= " -fastest" if $opt{'F'};
$opts .= $opt{'B'} ? " -boot $opt{B}" : " -nosupport";
$opts .= " -threads $opt{j}" if $opt{'e'} eq 'VeryFastTree';
msg("Opts: $opts");

run(
  "OMP_NUM_THREADS=$opt{j}"
 ." time $ftree $opts"
 ." -out $tmp.nwk1"
 ." -log $tmp.log"
 ." -nt $tmp.afa3"
);
copy("$tmp.log", "$out.$ftree.log");
copy("$tmp.nwk1", "$out.$ftree.nwk");

msg("Fixing branch lengths...");

run("gotree resolve < $tmp.nwk1 > $tmp.nwk2");
run(
  "raxml-ng --threads $cpus"
 ." --force perf_threads" # avoid bail out
 ." --evaluate --model GTR+G4 --blmin $BLMIN"
 ." --prefix $tmp"
 ." --tree $tmp.nwk2"
 ." --msa $tmp.afa3"
 ." --site-weights $tmp.afa3.weights"
);
#copy("$tmp.raxml.bestTreeCollapsed", "$tmp.nwk3");
copy("$tmp.raxml.bestTree", "$tmp.nwk3");

msg("Cleaning branches...");
run(
  " cat $tmp.nwk3"
 ." | gotree brlen round -p $PRECISION"
 ." | gotree collapse length -l 0"
 ." > $tmp.nwk4"
);

msg("Repopulating, sorting, rerooting tree...");

run(
  "gotree repopulate -i $tmp.nwk4 -g $tmp.dup"
 ." | nw_order -c a -"
 ." | nw_reroot - $outgroup"
 ." | gotree brlen scale -f 1"  
 ." > $tmp.nwk5"
);

  
msg("Collating output files");
#copy("$tmp.nwk1", "$out.raw.nwk");
copy("$tmp.nwk5", "$out.nwk");

run("nw_stats $out.nwk");

kovid_exit();

#......................................
sub num_lines {
  my($fname) = @_;
  my($L) = qx"wc -l < '$fname'";
  chomp $L;
  return $L;
}

#....................................................................
sub fasta2array {
  exe("seqkit");
  my $seq;
  for my $fname (@_) {
    msg("Loading FASTA: $fname");
    open my $TAB, '-|', "seqkit fx2tab --only-id '$fname'";
    while (<$TAB>) {
      chomp;
      push @$seq, [ split m/\t/ ];
      # replace empty seq with single N
      $seq->[-1] //= 'N';
    }
    close $TAB;
  }
  msg("Loaded", scalar(@$seq), "sequences.");
  return $seq;
}

#....................................................................
sub tmpdir {
  my($change) = @_;
  # CLEANUP->1 should be default for this
  my $dir = File::Temp->newdir();
  $dir or err("Can't create temp folder: $!");
  msg("Created tempdir: $dir");
  if ($change) {
    msg("Changing to tempdir: $dir");
    $old_dir = getcwd;
    chdir $dir;
  }
  return $dir;
}

#....................................................................
sub run { 
  msg("Running: @_"); 
  system(@_)==0 or err("Could not run: @_"); 
}