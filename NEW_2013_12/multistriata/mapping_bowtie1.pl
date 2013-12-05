#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

# VERSION: 0.2

##############
# PARAMETERS #
##############
my $folder = '/media/LOCAL_DATA_2/srobusta/splitted_originl_data/mating'; # complete path
my $index = '/media/LOCAL_DATA_2/srobusta/index/sr_genome_part_2013_10'; # complete path
#my $index = '/media/LOCAL_DATA_2/srobusta/index/sr_transcriptome_assembly_2013_10'; # complete path
my $outdir = '/media/LOCAL_DATA_2/srobusta/splitted_originl_data/mating_mapping_genome_out'; # complete path - IT MUST EXIST: WILL NOT BE CREATED
my $core = 12;
my $chunkmbs = 10240;
my $maxins = 500;
my $trim5 = 20;
my $trim3 = 20;
my $seedlen = 20;
my $nofw = 1;
my $norc = 0;

my $pm = new Parallel::ForkManager($core);
$pm->run_on_finish(
  sub {
    my($pid,$exit_code) = @_;
    print "** Just got out of the pool ".
          "with PID $pid and exit code: $exit_code\n";
  }
);

# Go into the directory containing the fastq files
chdir($folder) or die "\nCannot change in directory $folder: $!\n\n";

# Get all fastq files inside the folder
my @fastq = glob("*.fastq");

print "Found the following fastq files:\n\n";
print "$_\n" foreach @fastq;
print "\n\n";

# STEP 1: run bowtie
my $seen = {};
foreach my $fastq(@fastq) {
  my $new = "$fastq";
  $new =~ s/\_[12].fastq$//;
  next if exists $seen->{$new};
  $seen->{$new} ++;
  my $f1 = "$new\_1.fastq";
  my $f2 = "$new\_2.fastq";
  my $stranded;
  my $samout;
  my $log;

  if($nofw) {
    $samout = "$new\_minus\.sam";
    $log = "bowtie_$new\_minus\.log";
    $stranded = '--nofw';
  }
  elsif($norc) {
    $samout = "$new\_plus\.sam";
    $log = "bowtie_$new\_plus\.log";
    $stranded = '--norc';
  }
  else{
    $samout = "$new\.sam";
    $log = "bowtie_$new\.log";
    $stranded = ' ';
  }

  $samout = "$outdir\/$samout";

  my $command = "bowtie $index -t -q -1 $f1 -2 $f2 -p $core --chunkmbs $chunkmbs --maxins $maxins ".
                "--trim5 $trim5 --trim3 $trim3 --seedlen $seedlen --tryhard -a";

  $command .= " -S $samout $stranded > $log";

  run_command($command);
}

print "\n\nDONE BOTWIE STEP!!!\n\n\n";

# STEP 2: samtools
chdir($outdir) or die "\nCannot change in directory $folder: $!\n\n";
my @sam = glob("*.sam");
foreach my $sam(@sam) {
  # Forks and returns the pid for the child:
  my $pid = $pm->start and next;

  # Here is the parallelized block
  # -----------
  analyze_mapping($sam);
  # -----------

  # Terminates the child process
  $pm->finish;
}

sub analyze_mapping {

  my $sam = shift;
  my $name = "$sam";
  $name =~ s/\.sam$//;
  my $bam = "$name\.bam";
  my $sorted = "$name\.sorted";
  my $sortedbam = "$sorted\.bam";
  my $counts = "$name\.counts";

  my $command1 = "samtools view -b -S $sam > $bam 2> samtools_view_$name\.log";
  run_command($command1);
  run_command("rm $sam"); # to rescue space
  my $command2 = "samtools sort $bam $sorted > samtools_sort_$name\.log";  
  run_command($command2);
  run_command("rm $bam");# to rescue space
  my $command3 = "samtools index $sortedbam > samtools_index_$name\.log";
  run_command($command3);
  my $command4 = "samtools idxstats $sortedbam > $counts 2> samtools_idxstats_$name\.log";
  run_command($command4);

}

sub run_command {
  my $command = shift;
  print "\nLAUNCHING SYSTEM CALL:\n\t$command\n";
  system($command);
  die "ERROR using command:\n\t$command\:\n\t$!" unless $? == 0;
  print "$command DONE!\n";
}

__END__

# Example of a successful run from the command line

# RUN BOWTIE
# Every analysis can use more CPUs
bowtie index/pm_transcriptome_assembly_2013_10 -q -1 CIIO_1.fastq -2 CIIO_2.fastq -p 20 -S CIIO.sam --chunkmbs 10240 -t --maxins 500 --trim5 20 --trim3 20 --seedlen 20 --tryhard -a
bowtie index/pm_transcriptome_assembly_2013_10 -q -1 CIIP_1.fastq -2 CIIP_2.fastq -p 20 -S CIIP.sam --chunkmbs 10240 -t --maxins 500 --trim5 20 --trim3 20 --seedlen 20 --tryhard -a
bowtie index/pm_transcriptome_assembly_2013_10 -q -1 HATT_1.fastq -2 HATT_2.fastq -p 20 -S HATT.sam --chunkmbs 10240 -t --maxins 500 --trim5 20 --trim3 20 --seedlen 20 --tryhard -a
bowtie index/pm_transcriptome_assembly_2013_10 -q -1 HCUH_1.fastq -2 HCUH_2.fastq -p 20 -S HCUH.sam --chunkmbs 10240 -t --maxins 500 --trim5 20 --trim3 20 --seedlen 20 --tryhard -a
bowtie index/pm_transcriptome_assembly_2013_10 -q -1 HCUN_1.fastq -2 HCUN_2.fastq -p 20 -S HCUN.sam --chunkmbs 10240 -t --maxins 500 --trim5 20 --trim3 20 --seedlen 20 --tryhard -a
bowtie index/pm_transcriptome_assembly_2013_10 -q -1 HCUO_1.fastq -2 HCUO_2.fastq -p 20 -S HCUO.sam --chunkmbs 10240 -t --maxins 500 --trim5 20 --trim3 20 --seedlen 20 --tryhard -a

# CONVERT SAM TO BAM
# Every analysis use a single CPU
nohup samtools view -b -S CIIO.sam > CIIO.bam 2>nohup.ciio.bam &
nohup samtools view -b -S CIIP.sam > CIIP.bam 2>nohup.ciip.bam &
nohup samtools view -b -S HATT.sam > HATT.bam 2>nohup.hatt.bam &
nohup samtools view -b -S HCUH.sam > HCUH.bam 2>nohup.hcuh.bam &
nohup samtools view -b -S HCUN.sam > HCUN.bam 2>nohup.hcun.bam &
nohup samtools view -b -S HCUO.sam > HCUO.bam 2>nohup.hcuo.bam &

# SORT BAM
# Every analysis use a single CPU
nohup samtools sort CIIO.bam CIIO.sorted > nohup.ciio.sort &
nohup samtools sort CIIP.bam CIIP.sorted > nohup.ciip.sort &
nohup samtools sort HATT.bam HATT.sorted > nohup.hatt.sort &
nohup samtools sort HCUH.bam HCUH.sorted > nohup.hcuh.sort &
nohup samtools sort HCUN.bam HCUN.sorted > nohup.hcun.sort &
nohup samtools sort HCUO.bam HCUO.sorted > nohup.hcuo.sort &

# INDEX BAM
# Every analysis use a single CPU
nohup samtools index CIIO.sorted.bam > nohup.ciio.index &
nohup samtools index CIIP.sorted.bam > nohup.ciip.index &
nohup samtools index HATT.sorted.bam > nohup.hatt.index &
nohup samtools index HCUH.sorted.bam > nohup.hcuh.index &
nohup samtools index HCUN.sorted.bam > nohup.hcun.index &
nohup samtools index HCUO.sorted.bam > nohup.hcuo.index &

# COUNT MAPPINGS
# Every analysis use a single CPU
nohup samtools idxstats CIIO.sorted.bam > CIIO.counts 2>nohup.ciio.count &
nohup samtools idxstats CIIP.sorted.bam > CIIP.counts 2>nohup.ciip.count &
nohup samtools idxstats HATT.sorted.bam > HATT.counts 2>nohup.hatt.count &
nohup samtools idxstats HCUH.sorted.bam > HCUH.counts 2>nohup.hcuh.count &
nohup samtools idxstats HCUN.sorted.bam > HCUN.counts 2>nohup.hcun.count &
nohup samtools idxstats HCUO.sorted.bam > HCUO.counts 2>nohup.hcuo.count &

