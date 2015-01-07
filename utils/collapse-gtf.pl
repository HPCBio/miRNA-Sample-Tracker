#!/usr/bin/env perl
use 5.014;
use strict;
use warnings;
use autodie;
use Data::Dumper;
use Getopt::Long;
use Text::CSV_XS;
use lib $ENV{GIT_REPO}."/hpcbio/smalheiser_scripts/v2/lib";
use FeatureUtils;

my $features;
my $kc;
my $exon = 1;
my $source;

my $USAGE = <<USAGE;

    perl $0 -f <ANNOTATION.bed/gff> -s refseq

Requires -f

Reads a GTF file and writes (to STDOUT) a simplified GTF only giving start/end
downstream analyses.  Requires a specialized mapping of the knownCanonical
data to the kgXref tables for UCSC hg19.

USAGE

GetOptions(
    'features=s'    => \$features,
    'source=s'      => \$source,
);

open( my $fh, '<', $features) or die $!;

my $ct = 0;

my %genes2trx;

my %trx;

while (<$fh>) {
    chomp;
    my ($ref, $origsource, $type, $start, $end, $score, $strand, $phase, $attstr) = split("\t", $_);
    
    $source //= $origsource;
    
    next unless $type eq 'exon';
    
    my $atts = FeatureUtils::gtf_attributes($attstr);
    
    my $nm = $ref.':'.$atts->{gene_id}[0];
    
    # marshall this so genes with same names on different refs will still be okay
    push @{ $genes2trx{ $nm }{ $atts->{transcript_id}[0] } }, [$ref, $source, $type, $start, $end, $score, $strand, $phase, $atts];
}

for my $gene (sort keys %genes2trx) {
    my @trx;
    my @trx_names;
    my ($gene_start, $gene_end, $gene_strand);
    my ($realref, $realgene) = split(':', $gene, 2);
    for my $trx (sort keys %{$genes2trx{$gene}}) {
        my $trxdata = $genes2trx{$gene}{$trx};
        my $ref;
        my ($trx_ref, $trx_start, $trx_end, $trx_strand, $trx_atts);
        for my $exon (@{ $trxdata } ) {
            $trx_strand //= $exon->[6];
            $trx_ref //= $exon->[0];
            $trx_atts //= $exon->[8];
            if ($trx_ref ne $exon->[0]) {
                die "Trx with different references in exons, uh oh: $trx_ref != ".$exon->[0]." in gene $gene, transcript $trx";
            }
            $trx_start = $exon->[3] if !$trx_start || $exon->[3] < $trx_start;
            $trx_end = $exon->[4] if !$trx_end || $exon->[4] > $trx_end;
        }
        
        $gene_strand //= $trx_strand;
        push @trx, [$trx_ref, $trx_start, $trx_end, $trx_strand, $trx_atts];
        push @trx_names, $trx;
        $gene_start = $trx_start if !$gene_start || $trx_start < $gene_start;
        $gene_end = $trx_end if !$gene_end || $trx_end > $gene_end;
    }
    my $gene_atts = join('; ', 'gene_id "'.$realgene.'"', 'transcript_id "'.join(',', @trx_names).'"').';';
    say join("\t",
             $realref,
             $source,
             'gene',
             $gene_start,
             $gene_end,
             '.',
             $gene_strand,
             '.',
             $gene_atts
             );
    for my $trx (sort {$a->[1] <=> $b->[1]} @trx) {
        my $trx_atts = join('; ', map { $_.' "'.join(',', @{ $trx->[4]->{$_} } ) } sort keys %{ $trx->[4] }).';';

        say join("\t",
             $trx->[0],
             $source,
             'mRNA',
             $trx->[1],
             $trx->[2],
             '.',
             $trx->[3],
             '.',
             $trx_atts
             );
    }
}

