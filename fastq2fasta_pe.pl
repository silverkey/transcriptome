#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

# VERSION: 0.1

##################
#   PARAMETERS   #
##################
my $fastq_folder = '../../splitted_original_data/';
my $fastq1 = 'CIIP.1755.3.1605_1.fastq';
my $fastq2 = 'CIIP.1755.3.1605_2.fastq';

my $wd = cwd();
chdir($fastq_folder);

fastq2fasta_pe($fastq1,$fastq2);
print "fastq2fasta_pe.pl SUCCESS!\n";

sub fastq2fasta_pe {
  my $fastq1 = shift;
  my $fastq2 = shift;

  my $fasta1 = "$fastq1";
  $fasta1 =~ s/\.fastq$//;
  $fasta1 =~ s/\.fq$//;
  $fasta1 .= '.fasta';

  my $fasta2 = "$fastq2";
  $fasta2 =~ s/\.fastq$//;
  $fasta2 =~ s/\.fq$//;
  $fasta2 .= '.fasta';

  open(INL,$fastq1);
  open(INR,$fastq2);

  open(OUTL,">$fasta1");
  open(OUTR,">$fasta2");

  my $lrow;
  my $rrow;
  my $lhr = {};
  my $rhr = {};

  while($lrow = <INL>) {
    $rrow = <INR>;
    chomp($lrow);
    chomp($rrow);

    $lhr->{1} = $lrow;
    $lhr->{2} = <INL>;
    $lhr->{3} = <INL>;
    $lhr->{4} = <INL>;
  
    $rhr->{1} = $rrow;
    $rhr->{2} = <INR>;
    $rhr->{3} = <INR>;
    $rhr->{4} = <INR>;
  
    my $left_id = get_id($lrow);
    my $right_id = get_id($rrow);

    if($left_id ne $right_id) {
      print "\nERROR: Got a pair with unequal id: $left_id - $right_id\n";
    }

    $left_id .= '/1';
    my $fastaseq1 = to_fasta($left_id,$lhr->{2});
    print OUTL $fastaseq1;

    $right_id .= '/2';
    my $fastaseq2 = to_fasta($right_id,$rhr->{2});
    print OUTR $fastaseq2;

    $lrow = '';
    $rrow = '';
    $lhr = {};
    $rhr = {};
  }

  close(INL);
  close(INR);
  close(OUTL);
  close(OUTR);
}

sub get_id {
  my $row = shift;
  my @f = split(/\s+/,$row);
  my $id = $f[0];
  $id =~ s/\/[12]$//;
  return $id;
}

sub to_fasta {
  my ($id,$seq,$len) = @_;
  chomp($seq);
  # default to 80 characters of sequence per line
  $len = 80 unless $len;
  my $formatted_seq = ">$id\n";
  while (my $chunk = substr($seq, 0, $len, "")) {
    $formatted_seq .= "$chunk\n";
  }
  return $formatted_seq;
}
