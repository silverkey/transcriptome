#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

##############
# PARAMETERS #
##############
my $script = '/home/remo/src/trinityrnaseq_r2013-02-16/util/RSEM_util/run_RSEM_align_n_estimate.pl';
my $folder = '.';
my $transcriptome = '../../PM_transcriptome_no_norm.fasta';
my $seqtype = 'fa';
my $libtype = 'RF';
my $threads = 24;

my @seqio;
@seqio = glob('*.fastq') if $seqtype eq 'fq';
@seqio = glob('*.fasta') if $seqtype eq 'fa';

my $logfile = "$0";
$logfile =~ s/\.pl$/\.LOG/;
open(LOG,">$logfile");
print LOG "Found the following files:\n\n";
print LOG "$_\n" foreach @seqio;
print LOG "\n\n";

my $info = {};
foreach my $seqio(@seqio) {
  my $new = "$seqio";
  $new =~ s/\_[12]\.fastq$//;
  $new =~ s/\_[12]\.fasta$//;
  next if exists $info->{$new};
  $info->{$new}->{count} ++; 
  my $seqio1;
  my $seqio2;
  $seqio1 = "$new\_1.fastq" if $seqtype eq 'fq';
  $seqio1 = "$new\_1.fasta" if $seqtype eq 'fa';
  $seqio2 = "$new\_2.fastq" if $seqtype eq 'fq';
  $seqio2 = "$new\_2.fasta" if $seqtype eq 'fa';
  die "\nERROR: cannot find file $seqio1 from experiment $new\n" unless(-e $seqio1);
  die "\nERROR: cannot find file $seqio2 from experiment $new\n" unless(-e $seqio2);
  $info->{$new}->{left} = $seqio1;
  $info->{$new}->{right} = $seqio2;
  print LOG "Found experiment $new in files:\n$seqio1\n$seqio2\n\n";
}

print Dumper $info;

foreach my $exp(keys %$info) {
  my $left = $info->{$exp}->{left};
  my $right = $info->{$exp}->{right};

  my $cmd = "perl $script ";
  $cmd .= "--transcripts $transcriptome ";
  $cmd .= "--seqType $seqtype ";
  $cmd .= "--left $left ";
  $cmd .= "--right $right ";
  $cmd .= "--prefix $exp ";
  $cmd .= "--SS_lib_type $libtype "; 
  $cmd .= "--thread_count $threads";

  print LOG "\nLAUNCHING:\n$cmd\n";
  print "\nLAUNCHING:\n$cmd\n";
  my $ret = system($cmd);
  if($ret) {
    die "Error, cmd: $cmd died with ret: $ret";
  }
  print LOG "DONE!....\n";
}
