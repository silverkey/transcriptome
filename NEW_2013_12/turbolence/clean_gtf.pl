#!/usr/bin/perl
use strict;
use warnings;

my $usage = "\n\tUsage: perl $0 [gtf file]\n\n";
die $usage unless -e $ARGV[0];
my $in = $ARGV[0];
my $out = "$in";
$out =~ s/.gtf$/_pc_clean.gtf/;

open(IN,$in);
open(OUT,">$out");
while(my $line = <IN>) {
  chomp($line);
  my $row = "$line";
  $row =~ s/ $//g;
  $row =~ s/\;$//g;
  my @f = split(/\t/,$row);
  next unless $f[1] eq 'protein_coding';
  my @a = split(/\; /,$f[8]);
  foreach my $a(@a) {
    if($a =~ /\;/) {
      my $new = "$a";
      $new =~ s/\;/\-/g;
      print "$a will be changed in $new\n";
      print "Here is the original row:\n$line\n";
      $line =~ s/$a/$new/;
      print "Here is the new row:\n$line\n\n";
    }
  }
  print OUT "$line\n";
}
