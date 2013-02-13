#!/usr/bin/perl
use strict;
use warnings;

my $usage = "\nUSAGE: perl $0 [fastq file]\n\n";
die $usage unless scalar(@ARGV) == 1;
die $usage unless -e $ARGV[0];
die $usage unless $ARGV[0] =~ /.fastq$/;

my $fastq = $ARGV[0];

my $fastq1 = $fastq;
$fastq1 =~ s/\.fastq/\_1\.fastq/;

my $fastq2 = $fastq;
$fastq2 =~ s/\.fastq/\_2\.fastq/;

open(IN,$fastq);
open(OUT1,">$fastq1");
open(OUT2,">$fastq2");

my $row;
my $left = {};
my $right = {};

while($row = <IN>) {

  $left->{1} = $row;
  $left->{2} = <IN>;
  $left->{3} = <IN>;
  $left->{4} = <IN>;
  
  $right->{1} = <IN>;
  $right->{2} = <IN>;
  $right->{3} = <IN>;
  $right->{4} = <IN>;
  
  my $left_id = get_id($left->{1});
  my $right_id = get_id($right->{1});

  if($left_id ne $right_id) {
    die "\nERROR: Got a pair with unequal id: $left_id - $right_id\n";
  }

  print OUT1 $left->{1};
  print OUT1 $left->{2};
  print OUT1 $left->{3};
  print OUT1 $left->{4};

  print OUT2 $right->{1};
  print OUT2 $right->{2};
  print OUT2 $right->{3};
  print OUT2 $right->{4};

  $row = '';
  $left = {};
  $right = {};

}

sub get_id {
  my $row = shift;
  my @f = split(/\s+/,$row);
  my $id = $f[0];
  $id =~ s/\/[12]$//;
  return $id;
}
