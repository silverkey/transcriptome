#!/usr/bin/perl
use strict;
use warnings;

my $in = $ARGV[0];
open(IN,$in);

print "transcript\tchr\tstart\tend\tstrand\n";

while(my $row = <IN>) {
  chomp($row);
  my @f = split("\t",$row);
  next unless $f[2] eq 'mRNA';
  my $id = $f[8];
  $id =~ s/^ID\=(.+)\.\d+$/$1/;
  my $chr = $f[0];
  my $start = $f[3];
  my $end = $f[4];
  my $strand = $f[6];
  print "$id\t$chr\t$start\t$end\t$strand\n";
}

