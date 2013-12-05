#!/usr/bin/perl
use strict;
use warnings;

#my $command1 = '/home/remo/src/bwa-0.7.5a/bwa mem -t 20 genomes/allgenomes_april_2012_v6_one_per_species.fa fastq/34-1-1_1.fastq fastq/34-1-1_2.fastq > bwa_aln_34_pe.sam 2>ERR_bwa_aln_34_pe';
#my $command2 = '/home/remo/src/bwa-0.7.5a/bwa mem -t 20 genomes/allgenomes_april_2012_v6_one_per_species.fa fastq/35-1-1_1.fastq fastq/35-1-1_2.fastq > bwa_aln_35_pe.sam 2>ERR_bwa_aln_35_pe';

my $command1 = '/home/remo/src/bwa-0.7.5a/bwa mem -t 24 /home/remo/ANALYSIS/turbolence/genome/Phaeodactylum_tricornutum.ASM15095v2.20.dna.toplevel.fa /home/remo/ANALYSIS/turbolence/fastq/C2F70ACXX_sB2Aug13_13s006178-1-1_Ferrante_lane713s006178_sequence.txt > sB_bwa.sam 2>ERR_sB_bwa';
my $command2 = '/home/remo/src/bwa-0.7.5a/bwa mem -t 24 /home/remo/ANALYSIS/turbolence/genome/Phaeodactylum_tricornutum.ASM15095v2.20.dna.toplevel.fa /home/remo/ANALYSIS/turbolence/fastq/C2F70ACXX_sD2Aug13_13s006180-1-1_Ferrante_lane713s006180_sequence.txt > sD_bwa.sam 2>ERR_sD_bwa';
my $command3 = '/home/remo/src/bwa-0.7.5a/bwa mem -t 24 /home/remo/ANALYSIS/turbolence/genome/Phaeodactylum_tricornutum.ASM15095v2.20.dna.toplevel.fa /home/remo/ANALYSIS/turbolence/fastq/C2F70ACXX_sF2Aug13_13s006182-1-1_Ferrante_lane713s006182_sequence.txt > sF_bwa.sam 2>ERR_sF_bwa';
my $command4 = '/home/remo/src/bwa-0.7.5a/bwa mem -t 24 /home/remo/ANALYSIS/turbolence/genome/Phaeodactylum_tricornutum.ASM15095v2.20.dna.toplevel.fa /home/remo/ANALYSIS/turbolence/fastq/C2F70ACXX_tA1Aug13_13s006177-1-1_Ferrante_lane713s006177_sequence.txt > tA_bwa.sam 2>ERR_tA_bwa';
my $command5 = '/home/remo/src/bwa-0.7.5a/bwa mem -t 24 /home/remo/ANALYSIS/turbolence/genome/Phaeodactylum_tricornutum.ASM15095v2.20.dna.toplevel.fa /home/remo/ANALYSIS/turbolence/fastq/C2F70ACXX_tC2Aug13_13s006179-1-1_Ferrante_lane713s006179_sequence.txt > tC_bwa.sam 2>ERR_tC_bwa';
my $command6 = '/home/remo/src/bwa-0.7.5a/bwa mem -t 24 /home/remo/ANALYSIS/turbolence/genome/Phaeodactylum_tricornutum.ASM15095v2.20.dna.toplevel.fa /home/remo/ANALYSIS/turbolence/fastq/C2F70ACXX_tE2Aug13_13s006181-1-1_Ferrante_lane713s006181_sequence.txt > tE_bwa.sam 2>ERR_tE_bwa';

exec_command($command1);
exec_command($command2);
exec_command($command3);
exec_command($command4);
exec_command($command5);
exec_command($command6);

sub exec_command {
  my $command = shift;
  print "\nLAUNCHING SYSTEM CALL:\n\t$command\n";
  system($command);
  die "ERROR using command:\n\t$command\:\n\t$!" unless $? == 0;
  print "DONE!\n";
}

