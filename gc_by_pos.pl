#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use Data::Dumper;

sub countgc;
sub printpos;

my $bygene = 0; #1 if for each gene, 0 if for each genome
my $sequence;
my @position;
my %pos;
while (<>) {
  chomp;
  if ($_ =~ m/^>/) {
    if ($sequence) {
      if (length($sequence)%3==0) {
	countgc(\$sequence);
	printpos() if ($bygene);
      }
    }
    $sequence = '';
  } else {
    $sequence .= $_; 
  }
}

countgc(\$sequence);
printpos();

sub countgc {
  $sequence =~ tr/CT/GA/;
  for (my $i=0; $i < length($sequence); $i++) {
    $pos{$i%3}{substr $sequence, $i, 1}++;
  }
  return 1;
}

sub printpos {

  my @gc;
  @gc = map { $gc[$_] = $pos{$_}{'G'} } 0..2;
  my @at;
  @at = map { $at[$_] = $pos{$_}{'A'} } 0..2;
  
  print $gc[0]/($gc[0]+$at[0]),
  " ",$gc[1]/($gc[1]+$at[1]),
  " ",$gc[2]/($gc[2]+$at[2]),
  "\n";

  %pos = ();

  return 1;
}
