#!/usr/bin/perl
use strict;
use warnings;

# VERSION: 0.4

# The program collect all the files with .fastq extension in the folder and work on them.
# In case of PE already splitted files they MUST be called with the '_1' or '_2' before
# the .fastq extension.

##############
# PARAMETERS #
##############
my $folder = '.';
my $run_split_fastq = 0;
my $split_fasta = 'split_unique_pe_fastq.pl';
my $trimmomatic = '/home/remo/src/Trimmomatic-0.22/trimmomatic-0.22.jar';
my $illumina_adapters = 'illumina_adapters.fa';
my $trimmomatic_threads = '20';
my $trimmomatic_clip = '2:40:15';
my $trimmomatic_leading = '5';
my $trimmomatic_trailing = '5';
my $trimmomatic_sliding = '5:20';
my $trimmomatic_minlen = '30';
my $name_left = 'shrimp_filt_left';
my $name_right = 'shrimp_filt_right';
my $output_format = 'fasta';
my $seqfile_ending = '_sequence.txt';

# The program only accepts as output fasta or fastq
if($output_format eq 'fastq') {
  $name_left .= '.fastq';
  $name_right .= '.fastq';
}
elsif($output_format eq 'fasta') {
  $name_left .= '.fasta';
  $name_right .= '.fasta';
}
else {
  die "\nERROR: output format $output_format not recognized [fastq|fasta]\n";
}

# Hashref to populate with names and info on the fastq files
my $info = {};

# Go into the directory containing the fastq files
chdir($folder) or die "\nCannot change in directory $folder\n\n";

# Open LOG file
open(LOG,">preprocess_fastq_for_trinity.LOG");

# Get all fastq files inside the folder
my @fastq = glob("*$seqfile_ending");

print LOG "Found the following fastq files:\n\n";
print LOG "$_\n" foreach @fastq;
print LOG "\n\n";

# STEP 1: Split fastq in case PE are in the same file
my $seen = {};
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

    my $new = "$fastq";
    $new =~ s/$seqfile_ending$//;

    $info->{$new}->{sp_l} = $fastq1;
    $info->{$new}->{sp_r} = $fastq2;

    print LOG "Splitted $fastq ($new) in:\n$fastq1\n$fastq2\n\n";
  }
  else {
    my $new = "$fastq";
    $new =~ s/\_[12]$seqfile_ending$//;
    next if exists $seen->{$new};
    $seen->{$new} ++; 
    my $fastq1 = "$new\_1$seqfile_ending";
    my $fastq2 = "$new\_2$seqfile_ending";
    $info->{$new}->{sp_l} = $fastq1;
    $info->{$new}->{sp_r} = $fastq2;
    print LOG "Found lane $new in files:\n$fastq1\n$fastq2\n\n";
  }
}

# STEP 2: Run Trimmomatic
$seen = {};
foreach my $fastq(@fastq) {

  my $new = "$fastq";
  $new =~ s/\_[12]$seqfile_ending$//;
  next if exists $seen->{$new};
  $seen->{$new} ++; 

  my $l = $info->{$new}->{sp_l};
  my $r = $info->{$new}->{sp_r};

  # Check CASAVA version
  my $lcasava18 = check_casava18($l);
  my $rcasava18 = check_casava18($r);

  # Define names for trimmomatic outputs
  my $lp = $l;
  $lp =~ s/$seqfile_ending/\_trim_paired.fastq/;
  my $lunp = $l;
  $lunp =~ s/$seqfile_ending/\_trim_unpaired.fastq/;
  my $rp = $r;
  $rp =~ s/$seqfile_ending/\_trim_paired.fastq/;
  my $runp = $r;
  $runp =~ s/$seqfile_ending/\_trim_unpaired.fastq/;

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

# STEP 3: Prepare for Digital Normalization and/or for Trinity however converting the fastq data in fasta
# According to the Trinity website/manual we must be sure that in the case of paired reads, the left read
# names end with suffix /1 and the right read names end with /2.

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

    if($output_format eq 'fastq') {
      print LEFT "$left_id\/1\n";
      print LEFT $lhr->{2};
      print LEFT $lhr->{3};
      print LEFT $lhr->{4};

      print RIGHT "$right_id\/2\n";
      print RIGHT $rhr->{2};
      print RIGHT $rhr->{3};
      print RIGHT $rhr->{4};
    }

    elsif($output_format eq 'fasta') {
      $left_id .= '/1';
      my $fastaseq1 = to_fasta($left_id,$lhr->{2});
      print LEFT $fastaseq1;

      $right_id .= '/2';
      my $fastaseq2 = to_fasta($right_id,$rhr->{2});
      print RIGHT $fastaseq2;
    }

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
  if($output_format eq 'fasta') {
    $id = s/^\@//;
  }
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
