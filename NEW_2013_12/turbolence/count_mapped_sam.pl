#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

#--------------------------------------------------------------------------
# THIS SCRIPT DO THE FOLLOWING:
# 1) ENTER INTO THE FOLDER WITH THE SAM TO ANALYZE
# 2) RUN HTSEQ WITH THE GIVEN PARAMETER AGAINST THE GIVEN GFF
#--------------------------------------------------------------------------

#------------------------------------------
# NUMBER OF THREADS
#------------------------------------------
my $threads = 24;

#------------------------------------------
# PARAMETERS ASSOCIATE TO FILES AND FOLDERS
#------------------------------------------
my $samdir = '.'; # the folder containing the bam
my $gff = '/home/remo/ANALYSIS/turbolence/map/Phaeodactylum_tricornutum.ASM15095v2.20_pc_clean.gtf';
#------------------------------------------

#------------------------------------------
# PARAMETERS ASSOCIATED TO HTSEQ ANALYSES
#------------------------------------------
my $mode = 'union';
my $stranded = 'no';
my $type = 'exon';
my $idattr = 'gene_id';
#------------------------------------------

#------------------------------------------
# LET'S THE ANALYSES BEGIN...
#------------------------------------------
chdir($samdir);

my @sam = glob('*.sam');

# NOW WE RUN HTSEQ ON EACH OF THE BAM FILES IN PARALLEL
my $pm = new Parallel::ForkManager($threads);
$pm->run_on_finish(
  sub {
    my($pid,$exit_code) = @_;
    print "** Just got out of the pool ".
          "with PID $pid and exit code: $exit_code\n";
  }
);
foreach my $sam(@sam) {
  my $pid = $pm->start and next;
  my $htseq = "htseq-count --mode=$mode --stranded=$stranded --type=$type --idattr=$idattr $sam $gff > COUNTS_$sam 2>ERR_count_$sam";
  exec_command($htseq);
  $pm->finish;
}

# EXEC SYSTEM CALL AND SAFELY CHECK FOR ERRORS
sub exec_command {
  my $command = shift;
  print "\nLAUNCHING SYSTEM CALL:\n\t$command\n";
  system($command);
  die "ERROR using command:\n\t$command\:\n\t$!" unless $? == 0;
  print "DONE!\n";
}
