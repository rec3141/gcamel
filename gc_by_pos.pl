#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use Data::Dumper;

sub countgc;
sub printpos;
sub resample;

my $debug = 0;
my $bygene = shift(@ARGV); #0 if for whole file, 1 if for each gene
my $reps = shift(@ARGV); #number of bootstrap replicates to output
$bygene = 0 if ($reps > 1); # if bootstraps, only whole file is supported
my $header = $ARGV[0];
my $sequence;
my @position;
my %pos;
LINE: while (<>) {
  chomp;
  #start a new gene
  if ($_ =~ m/^>/) {
    #is there a sequence already?
    #if so, it needs processing
    if ($sequence) {
      #is it the right length?
      #if not, skip it
      if (length($sequence)%3==0) {
	# wait until end to process if multiple reps
	next LINE if ($reps > 1);
	# otherwise run stats
	countgc(\$sequence);
	# but only print if need each gene processed
	printpos($header) if ($bygene == 1);
	}
      }
    $sequence = '' unless ($reps > 1);
    $header=$_ if ($bygene == 1);
  } else {
    $sequence .= $_; 
  }
}

if ($reps > 1) {
  resample($reps,\$sequence);
} else {
  countgc(\$sequence);
  printpos($header);
}

exit 0;


# SUBROUTINES

sub countgc {
  my $tmpseq = ${(shift)};
  my $seq = $tmpseq;
  $seq =~ s/[^ATGC]/0/g;
  $seq =~ tr/ATGC/0011/;
  my $invseq = $seq;
  $invseq =~ tr/10/01/;
  my $codons = length($seq)/3;
  my $mask1 = '100' x $codons;
  my $mask2 = '010' x $codons;
  my $mask3 = '001' x $codons;
  $pos{0}{'GC'} += (($mask1 & $seq) =~ tr/1//);
  $pos{1}{'GC'} += (($mask2 & $seq) =~ tr/1//);
  $pos{2}{'GC'} += (($mask3 & $seq) =~ tr/1//);
  $pos{0}{'AT'} += (($mask1 & $invseq) =~ tr/1//);
  $pos{1}{'AT'} += (($mask2 & $invseq) =~ tr/1//);
  $pos{2}{'AT'} += (($mask3 & $invseq) =~ tr/1//);
  
  return 1;
}

sub printpos {
  my $head = shift;
  my @gc;
  @gc = map { $gc[$_] = $pos{$_}{'GC'} } 0..2;
  my @at;
  @at = map { $at[$_] = $pos{$_}{'AT'} } 0..2;

  $debug>0 ? print "$head\t" : undef;

print $gc[0]/($gc[0]+$at[0]),
  "\t",$gc[1]/($gc[1]+$at[1]),
  "\t",$gc[2]/($gc[2]+$at[2]),
  "\n";
  
  %pos = ();

  return 1;
}

sub resample {
  my $times = shift;
  my $seq =  ${(shift)};

  my @codons = unpack("(A3)*", $seq);
  for (my $i=1; $i<=$times; $i++) {
    my $bootseq = '';
    for (my $j=1; $j < scalar(@codons); $j++) {
        $bootseq .= $codons[ rand(@codons) ];
    }
    countgc(\$bootseq);
    printpos($i . " " . $header);
  }
  
  return 1;
}