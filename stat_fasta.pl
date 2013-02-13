#!/usr/bin/perl
use strict;
use warnings;
use Statistics::Descriptive;
use Bio::SeqIO;
use Data::Dumper;

my $fasta = $ARGV[0];

my $usage = "\n\tUsage: perl $0 [fasta]\n\n";
die $usage unless scalar(@ARGV) == 1;
die $usage unless -e $fasta;

my $in = Bio::SeqIO->new(-file => $fasta,
                         -format => 'fasta');

my @length;
my @longest;
my @shortest;

my $n   = 0;
my $k01 = 0;
my $k02 = 0;
my $k05 = 0;
my $k1  = 0;
my $k2  = 0;
my $k3  = 0;
my $k5  = 0;
my $k10 = 0;

my $stat = Statistics::Descriptive::Full->new();

while(my $seq = $in->next_seq) {
  my $length = $seq->length;
  push(@length,$length);
  if($length > 10000) {
    push(@longest,$seq->id.' ---> '.$length);
    $k01 ++; $k02 ++; $k05 ++; $k1 ++; $k2 ++; $k3 ++; $k5 ++; $k10 ++;
  }
  elsif($length > 5000) {
    $k01 ++; $k02 ++; $k05 ++; $k1 ++; $k2 ++; $k3 ++; $k5 ++;
  }
  elsif($length > 3000) {
    $k01 ++; $k02 ++; $k05 ++; $k1 ++; $k2 ++; $k3 ++;
  }
  elsif($length > 2000) {
    $k01 ++; $k02 ++; $k05 ++; $k1 ++; $k2 ++;
  }
  elsif($length > 1000) {
    $k01 ++; $k02 ++; $k05 ++; $k1 ++;
  }
  elsif($length > 500) {
    $k01 ++; $k02 ++; $k05 ++; 
  }
  elsif($length > 200) {
    $k01 ++; $k02 ++;
  }
  elsif($length > 100) {
    push(@shortest,$seq->id.' ---> '.$length);
    $k01 ++;
  }
  $n ++;
}

$stat->add_data(@length);
my $mean = $stat->mean();
$mean =~ s/(\d+)\.\d+/$1/;
my $median = $stat->median();
$median =~ s/(\d+)\.\d+/$1/;
my $max = $stat->max();
my $min = $stat->min();

print "\nShortest transcripts are:\n";
print join("\n",@shortest);
print "\n\nLongest transcripts are:\n";
print join("\n",@longest);
print "\n\nTotal Number of Sequences = $n\n";
print "Average length = $mean\n";
print "Median length = $median\n";
print "Minimum length = $min\n";
print "Maximum length = $max\n";
print "Number of Sequences longer than 10000 = $k10\n";
print "Number of Sequences longer than 5000  = $k5\n";
print "Number of Sequences longer than 3000  = $k3\n";
print "Number of Sequences longer than 2000  = $k2\n";
print "Number of Sequences longer than 1000  = $k1\n";
print "Number of Sequences longer than 500   = $k05\n";
print "Number of Sequences longer than 200   = $k02\n";
print "Number of Sequences longer than 100   = $k01\n\n";
