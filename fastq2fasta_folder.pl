#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

# VERSION: 0.1

##################
#   PARAMETERS   #
##################
my $fastq_folder = '../../splitted_original_data/';

my $wd = cwd();
chdir($fastq_folder);

# Get all fastq files inside the folder
my @fastq = glob('*.fastq');

foreach my $fastq(@fastq) {
  print "Converting $fastq to fasta....\n";
  fastq2fasta_pe($fastq);
}
print "\nfastq2fasta_folder.pl SUCCESS!\n";

sub fastq2fasta_pe {
  my $fastq = shift;

  my $fasta = "$fastq";
  $fasta =~ s/\.fastq$//;
  $fasta =~ s/\.fq$//;
  $fasta .= '.fasta';

  open(IN,$fastq);
  open(OUT,">$fasta");

  my $row;
  my $hr = {};

  while($row = <IN>) {

    chomp($row);
    $hr->{1} = $row;
    $hr->{2} = <IN>;
    $hr->{3} = <IN>;
    $hr->{4} = <IN>;
  
    my $id = $row;

    my $fastaseq = to_fasta($id,$hr->{2});
    print OUT $fastaseq;

    $row = '';
    $hr = {};
  }

  close(IN);
  close(OUT);
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
