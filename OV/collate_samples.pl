#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

# IZ10_1_1.fastq
# IZ10_2_1.fastq
# IZ10_1_2.fastq
# IZ10_2_2.fastq

my @glob = glob('*.fastq');

my $href = {};

foreach my $file(@glob) {
  my $name = "$file";
  $name =~ s/\.fastq$//;
  my @info = split(/\_/,$name);
  $href->{$info[0]}->{$info[1]}->{$info[2]} = $file;
}

print Dumper $href;

foreach my $samp(keys %$href) {

  my $left1 = $href->{$samp}->{1}->{1};
  my $left2 = $href->{$samp}->{2}->{1};
  my $right1 = $href->{$samp}->{1}->{2};
  my $right2 = $href->{$samp}->{2}->{2};

  my $command1 = "cat $left1 \> $samp\_1.fastq";
  my $command2 = "cat $left2 \>> $samp\_1.fastq";
  my $command3 = "cat $right1 \> $samp\_2.fastq";
  my $command4 = "cat $right2 \>> $samp\_2.fastq";

  print "$command1\n";
  print "$command2\n";
  print "$command3\n";
  print "$command4\n";
  print "\n\n";
}
