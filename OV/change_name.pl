#!/usr/bin/perl
use strict;
use warnings;

# _IZ10_13s002471-1-1_Zarrella_

print "change_name.pl modifications:\n\n";

my @glob1 = glob('*_1_sequence.txt');

foreach my $file (@glob1) {
  my $new = "$file";
  $new =~ /^.+\_(IZ\d+)\_.+\-(\d)\-\d.+(\_1\_sequence\.txt)$/;
  $new = "$1\_$2$3";
  $new =~ s/\_sequence\.txt$/\.fastq/;
  print "mv $file $new\n";
  system("mv $file $new");
}

my @glob2 = glob('*_2_sequence.txt');

foreach my $file (@glob2) {
  my $new = "$file";
  $new =~ /^.+\_(IZ\d+)\_.+\-(\d)\-\d.+(\_2\_sequence\.txt)$/;
  $new = "$1\_$2$3";
  $new =~ s/\_sequence\.txt$/\.fastq/;
  print "mv $file $new\n";
  system("mv $file $new");
}

print "\n\n\n\n";
