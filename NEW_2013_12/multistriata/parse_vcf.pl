#!/usr/bin/perl
use strict;
use warnings;
use Vcf;
use Data::Dumper;

my $vcf = Vcf->new(file=>'transcriptome_mapping_SNP_calls.vcf');
$vcf->parse_header();

my @samples = $vcf->get_samples;
print_header(\@samples);

while (my $x = $vcf->next_data_hash()) {
  if($x->{INFO}->{DP} >= 100) {
    print_line($x);
    print_anno($x,\@samples);
  }
}

sub print_header {
  my $samples = shift;
  print "seqid\tpos\tdp\tqual\tref\talt";
  for(my $c=0;$c<scalar(@$samples);$c++) {
    my $sample = $samples[$c];
    $sample =~ s/.bam$//;
    print "\t$sample\_gt\t$sample\_dp\t$sample\_gq\t$sample\_pl";
  }
  print "\n";
}

sub print_line {
  my $x = shift;
  print $x->{CHROM}."\t";
  print $x->{POS}."\t";
  print $x->{INFO}->{DP}."\t";
  print $x->{QUAL}."\t";
  print $x->{REF}."\t";
  print join(',',@{$x->{ALT}});
}

sub print_anno {
  my $x = shift;
  my $samples = shift;
  for(my $c=0;$c<scalar(@$samples);$c++) {
    my $sample = $samples[$c];
    print "\t".$x->{gtypes}->{$sample}->{GT};
    print "\t".$x->{gtypes}->{$sample}->{DP};
    print "\t".$x->{gtypes}->{$sample}->{GQ};
    print "\t".$x->{gtypes}->{$sample}->{PL};
  }
  print "\n";
}

__END__

sample	mating	size	other
CIIP	-	S	nonna
HATT	-	L	nipote1
HCUN	-	S	nipote1
CIIO	+	S	nonno
HCUO	+	L	nipote2 (clone genoma)
HCUH	+	S	nipote2 (clone genoma)


$VAR1 = {
          'FORMAT' => [
                        'GT',
                        'PL',
                        'DP',
                        'GQ'
                      ],
          'QUAL' => '38.7',
          'ID' => '.',
          'CHROM' => 'comp25_c0_seq1',
          'INFO' => {
                      'DP' => '19',
                      'MQ' => '20',
                      'DP4' => '1,11,1,4',
                      'VDB' => '9.536029e-02',
                      'AF1' => '0.324',
                      'AC1' => '4',
                      'PV4' => '0.51,1,1,0.32',
                      'RPB' => '-4.743417e-01',
                      'FQ' => '40.3'
                    },
          'FILTER' => [
                        '.'
                      ],
          'gtypes' => {
                        'HCUN.bam' => {
                                        'GQ' => '3',
                                        'DP' => '0',
                                        'PL' => '0,0,0',
                                        'GT' => '0/0'
                                      },
                        'HATT.bam' => {
                                        'GQ' => '25',
                                        'DP' => '8',
                                        'PL' => '0,24,83',
                                        'GT' => '0/0'
                                      },
                        'HCUH.bam' => {
                                        'GQ' => '3',
                                        'DP' => '0',
                                        'PL' => '0,0,0',
                                        'GT' => '0/0'
                                      },
                        'HCUO.bam' => {
                                        'GQ' => '3',
                                        'DP' => '0',
                                        'PL' => '0,0,0',
                                        'GT' => '0/0'
                                      },
                        'CIIP.bam' => {
                                        'GQ' => '13',
                                        'DP' => '4',
                                        'PL' => '0,12,70',
                                        'GT' => '0/0'
                                      },
                        'CIIO.bam' => {
                                        'GQ' => '9',
                                        'DP' => '5',
                                        'PL' => '82,15,0',
                                        'GT' => '1/1'
                                      }
                      },
          'REF' => 'A',
          'ALT' => [
                     'T'
                   ],
          'POS' => '152'
        };

