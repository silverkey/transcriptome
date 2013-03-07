#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;

# VERSION: 2.0

my $fasta = $ARGV[0];
my $length = $ARGV[1];

my $usage = "\n\tUsage: perl $0 [fasta] [length]\n\n";
die $usage unless scalar(@ARGV) == 2;
die $usage unless -e $fasta;
die $usage unless $length >= 1;

my $new = "$fasta";
$new =~ s/\.(.+)$//;
$new .= "_comp_length_$length\.fasta";

my $nopass = "$new";
$nopass =~ s/\.fasta$/_NOPASS.fasta/;

my $href = {};
my $lengthok = {};
my $written = {};

my $in = Bio::SeqIO->new(-file => $fasta,
                         -format => 'fasta');

my $out = Bio::SeqIO->new(-file => ">$new",
                          -format => 'fasta');

my $ref = Bio::SeqIO->new(-file => ">$nopass",
                          -format => 'fasta');

while(my $seq = $in->next_seq) {
  next unless $seq->length >= $length;
  my $length = $seq->length;
  my $id = $seq->id;
  $lengthok->{$id} = $seq;
  my(@val) = split(/\_/,$id);
  my $compc = join('_',($val[0],$val[1]));
  my $at = $val[2];
  if(exists $href->{$compc}) {
    my $seen = $href->{$compc};
    if($seq->length > $seen->length) {
      $href->{$compc} = $seq;
    }
  }
  else {
    $href->{$compc} = $seq;
  }
}

foreach my $compc(keys %$href) {
  $written->{$href->{$compc}->id} ++;
  $out->write_seq($href->{$compc});
}

foreach my $id(keys %$lengthok) {
  next if exists $written->{$id};
  $ref->write_seq($lengthok->{$id});
}
