#!/usr/bin/perl
use strict;
use warnings;

# Parameters
my $folder = '.';
my $run_split_fastq = 1;
my $split_fasta = 'split_unique_pe_fastq.pl';
my $trimmomatic = '/home/remo/src/Trimmomatic-0.22/trimmomatic-0.22.jar';
my $illumina_adapters = 'illumina_adapters.fa';
my $trimmomatic_threads = '20';
my $trimmomatic_clip = '2:40:15';
my $trimmomatic_leading = '5';
my $trimmomatic_trailing = '5';
my $trimmomatic_sliding = '5:20';
my $trimmomatic_minlen = '100';
my $name_left = 'global_left.fastq';
my $name_right = 'global_right.fastq';

# Hashref to populate with names and info on the fastq files
my $info = {};

# Go into the directory containing the fastq files
chdir($folder) or die "\nCannot change in directory $folder\n\n";

# Open LOG file
open(LOG,">preprocess_fastq_for_trinity.LOG");

# Get all fastq files inside the folder
my @fastq = glob('*.fastq');

print LOG "Found the following fastq files:\n\n";
print LOG "$_\n" foreach @fastq;
print LOG "\n\n";

# STEP 1: Split fastq
foreach my $fastq(@fastq) {
  if($run_split_fastq) {
    print LOG "Launching perl $split_fasta $fastq\n";
    system("perl $split_fasta $fastq");

    # Do not change this as it is reflecting naming created
    # by the script splitting fastq unique pe files!
    my $fastq1 = $fastq;
    $fastq1 =~ s/\.fastq/\_1\.fastq/;
    my $fastq2 = $fastq;
    $fastq2 =~ s/\.fastq/\_2\.fastq/;

    $info->{$fastq}->{sp_l} = $fastq1;
    $info->{$fastq}->{sp_r} = $fastq2;

    print LOG "Splitted $fastq in:\n$fastq1\n$fastq2\n\n";
  }
}

# STEP 2: Run Trimmomatic
foreach my $fastq(@fastq) {

  my $l = $info->{$fastq}->{sp_l};
  my $r = $info->{$fastq}->{sp_r};

  # Check CASAVA version
  my $lcasava18 = check_casava18($l);
  my $rcasava18 = check_casava18($r);

  # Define names for trimmomatic outputs
  my $lp = $l;
  $lp =~ s/\.fastq/\_trim_paired.fastq/;
  my $lunp = $l;
  $lunp =~ s/\.fastq/\_trim_unpaired.fastq/;
  my $rp = $r;
  $rp =~ s/\.fastq/\_trim_paired.fastq/;
  my $runp = $r;
  $runp =~ s/\.fastq/\_trim_unpaired.fastq/;

  if($lcasava18 ne $rcasava18) {
    print LOG "ERROR: splitted fastq are not associated to a similar CASAVA version\n";
  }

  $info->{$fastq}->{casava18} = "$lcasava18\/$rcasava18";
  my $phred = 64;
  $phred = 33 if $lcasava18 eq 'T';

  my $trimmomatic_command = "java -classpath $trimmomatic org.usadellab.trimmomatic.TrimmomaticPE ";
  $trimmomatic_command .= "-threads $trimmomatic_threads -phred$phred ";
  $trimmomatic_command .= "$l $r $lp $lunp $rp $runp ";
  $trimmomatic_command .= "ILLUMINACLIP:$illumina_adapters\:$trimmomatic_clip ";
  $trimmomatic_command .= "LEADING:$trimmomatic_leading TRAILING:$trimmomatic_trailing SLIDINGWINDOW:$trimmomatic_sliding ";
  $trimmomatic_command .= "MINLEN:$trimmomatic_minlen";

  print LOG "Launching Trimmomatic with the following command:\n$trimmomatic_command\n\n";

  system($trimmomatic_command);

  $info->{$fastq}->{tr_l} = $lp;
  $info->{$fastq}->{tr_r} = $rp;
}

# STEP 3: Prepare for Digital Normalization
# According to the Trinity website/manual:
# Before running the normalization, be sure that in the case of paired reads, the left read names end with suffix /1 and the right read names end with /2.
open(LEFT,">$name_left");
open(RIGHT,">$name_right");

foreach my $fastq(@fastq) {

  my $lp = $info->{$fastq}->{tr_l};
  my $rp = $info->{$fastq}->{tr_r};

  open(INL,$lp);
  open(INR,$rp);

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
      print LOG "\nERROR: Got a pair with unequal id: $left_id - $right_id\n";
    }

    print LEFT "$left_id\/1\n";
    print LEFT $lhr->{2};
    print LEFT $lhr->{3};
    print LEFT $lhr->{4};

    print RIGHT "$right_id\/2\n";
    print RIGHT $rhr->{2};
    print RIGHT $rhr->{3};
    print RIGHT $rhr->{4};

    $lrow = '';
    $rrow = '';
    $lhr = {};
    $rhr = {};
  }

  close(INL);
  close(INR);
}

sub get_id {
  my $row = shift;
  my @f = split(/\s+/,$row);
  my $id = $f[0];
  $id =~ s/\/[12]$//;
  return $id;
}

sub check_casava18 {
  my $f = shift;
  my $casava18;
  open(IN,$f);
  my $row = <IN>;
  chomp($row);
  my @f = split(/\s+/,$row);
  my $el = scalar(@f);
  # CASAVA < 1.8 only single element or the second element without ':'
  if($el == 1) {
    print LOG "$row recognized as CASAVA < 1.8\n";
    $casava18 = 'F';
  }
  elsif($f[1] !~ /\:/) {
    print LOG "$row recognized as CASAVA < 1.8\n";
    $casava18 = 'F';
  }
  elsif($f[1] =~ /\:/) {
    print LOG "$row recognized as CASAVA >= 1.8\n";
    $casava18 = 'T';
  }
  else {
    print LOG "ERROR: Cannot recognize CASAVA format\n";
    $casava18 = 'NA';
  }
  close(IN);
  return $casava18;
}
