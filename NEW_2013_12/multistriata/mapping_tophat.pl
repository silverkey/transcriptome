#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

# VERSION: 0.2

##############
# PARAMETERS #
##############
#my $folder = '/home/remo/ANALYSIS/pmultistriata/FINAL/fastq'; # complete path
#my $index = '/home/remo/ANALYSIS/pmultistriata/FINAL/pmultistriata/genome/pm_genome_allpaths_2'; # complete path
#my $outdir = '/media/LOCAL_DATA_2/multistriata_genome_mapping'; # complete path - IT MUST EXIST: WILL NOT BE CREATED
my $folder = '/media/LOCAL_DATA_2/srobusta/splitted_originl_data/mating'; # complete path
my $index = '/media/LOCAL_DATA_2/srobusta/index/sr_genome_part_2013_10_B2'; # complete path
my $outdir = '/media/LOCAL_DATA_2/srobusta/splitted_originl_data/mating_mapping_genome_tophat_out'; # complete path - IT MUST EXIST: WILL NOT BE CREATED
my $core = 24;
my $library_type = 'fr-firststrand';

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
  my $log;

  my $out = "$outdir\/$new";

  my $command = "tophat --library-type $library_type --num-threads $core --output-dir $out $index $f1\,$f2";
  run_command($command);
}

print "\n\nDONE TOPHAT STEP!!!\n\n\n";

sub run_command {
  my $command = shift;
  print "\nLAUNCHING SYSTEM CALL:\n\t$command\n";
  system($command);
  die "ERROR using command:\n\t$command\:\n\t$!" unless $? == 0;
  print "$command DONE!\n";
}

__END__
