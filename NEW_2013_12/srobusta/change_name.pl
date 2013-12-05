#!/usr/bin/perl
use strict;
use warnings;

my @glob1 = glob('*_1.fasta');

foreach my $file (@glob1) {
  my $new = "$file";
  $new =~ /(\w+)\..+(\_1.fasta)/;
  $new = $1.$2;
  print "mv $file $new\n";
  system("mv $file $new");
}

my @glob2 = glob('*_2.fasta');

foreach my $file (@glob2) {
  my $new = "$file";
  $new =~ /(\w+)\..+(\_2.fasta)/;
  $new = $1.$2;
  print "mv $file $new\n";
  system("mv $file $new");
}
