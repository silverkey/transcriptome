#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;

my $transcriptome = 'PM_transcriptome_no_norm_comp_length_200.fasta';
my $results = 'genes.counts.matrix.minus_vs_plus.edgeR.DE_results';
my $counts = 'genes.counts.matrix';
my $analysis = 'mating';
my $type = 'genes';
my $minfdr = 0.05;

test_file($transcriptome);
test_file($results);
test_file($counts);

my $seqhr = populate_seq();
my($selhr,$selcol) = Rtab_with_rowname_and_head_to_href($results);
#print Dumper $selhr;
my($couhr,$coucol) = Rtab_with_rowname_and_head_to_href($counts);
#print Dumper $couhr;

open(OUT,">selected_$analysis\_$type\.xls");
print OUT "id\t";
foreach my $h(sort @$coucol) {
  print OUT "$h\t";
}
foreach my $h(sort @$selcol) {
  print OUT "$h\t";
}
print OUT "seq\n";

my $seqio = Bio::SeqIO->new(-file => ">selected_$analysis\_$type\.fa",
                            -format => 'fasta');

foreach my $id(keys %$selhr) {
  my $stat = $selhr->{$id};
  next unless $stat->{FDR} <= $minfdr;
  my $count = $couhr->{$id};
  my $seq = $seqhr->{$id};
  print OUT "$id\t";
  foreach my $sample(sort keys %$count) {
    print OUT $count->{$sample}."\t";
  }
  foreach my $val(sort keys %$stat) {
    print OUT $stat->{$val}."\t";
  }
  print OUT $seq->seq."\n";
  $seq->desc('');
  $seqio->write_seq($seq);
}

sub populate_seq {
  my $href = {};
  my $seqio = Bio::SeqIO->new(-file => $transcriptome,
                              -format => 'fasta');

  while(my $seq = $seqio->next_seq) {
    my $id = $seq->id;
    $id =~ s/\_seq\d+$// if $type eq 'genes';
    $href->{$id} = $seq;
  }
  return $href;
}

sub test_file {
  my $file = shift;
  print "File: $file\tOK\n" if -e $file;
  die "\nERROR: Cannot find file $file\n\n" unless -e $file;
}

sub Rtab_with_rowname_and_head_to_href {
  my $table = shift;
  my $href = {};
  open(IN,$table);
  my $head = <IN>;
  chomp $head;
  my @col = split(/\t/,$head);
  while(my $row = <IN>) {
    chomp $row;
    my @val = split(/\t/,$row);
    my $id = $val[0];
    for(my $x=1;$x<=$#val;$x++) {
      $href->{$id}->{$col[$x-1]} = $val[$x];
    }
  }
  return($href,\@col);
}
