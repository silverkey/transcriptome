#/usr/bin/perl
use strict;
use warnings;

# nohup perl /home/remo/src/trinityrnaseq_r2013-02-25/util/RSEM_util/merge_RSEM_frag_counts_single_table.pl YC1.isoforms.results YC2.isoforms.results YC3.isoforms.results YC4.isoforms.results YC5.isoforms.results YC6.isoforms.results YC7.isoforms.results YC8.isoforms.results > transcripts.counts.matrix &
# nohup perl /home/remo/src/trinityrnaseq_r2013-02-25/util/RSEM_util/merge_RSEM_frag_counts_single_table.pl YC1.genes.results YC2.genes.results YC3.genes.results YC4.genes.results YC5.genes.results YC6.genes.results YC7.genes.results YC8.genes.results > genes.counts.matrix &

my $trinity_dir = '/home/remo/src/trinityrnaseq_r2013-02-25/';
my $rsem_util_dir = 'util/RSEM_util/';
my $merge_script = 'merge_RSEM_frag_counts_single_table.pl';

my @isoform = glob('*isoforms.results');
my $isoform_command = $trinity_dir.$rsem_util_dir.$merge_script.' ';
foreach my $isoform(@isoform) {
  $isoform_command .= "$isoform ";
}
$isoform_command .= '> transcripts.counts.matrix';
print "\n$isoform_command\n\n";
system("$isoform_command");

my @gene = glob('*genes.results');
my $gene_command = $trinity_dir.$rsem_util_dir.$merge_script.' ' ;
foreach my $gene(@gene) {
  $gene_command .= "$gene ";
}
$gene_command .= '> genes.counts.matrix';
print "\n$gene_command\n\n";
system("$gene_command");


