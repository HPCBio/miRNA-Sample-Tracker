package FeatureUtils;

use 5.016;
use strict;
use warnings;
no if $] >= 5.017011, warnings => 'experimental::smartmatch';
use feature "switch";
use parent qw(Exporter);
use Text::CSV_XS;
our @EXPORT_OK=qw();

sub gtf_attributes {
    my ($data) = @_;
    my %atts = map {
        my ($key, $rest) = split(/\s/, $_, 2);
        $rest =~ s/['"]//g;
        $key => [split ',', $rest]
        } split(/\s*;\s*/, $data);
    \%atts;
}

# Reads UCSC kgXref data and creates a hash ref based on the key given (see below);
# if a specific ID isn't given, we assume for simplicity 'mRNA', as this
# corresponds to the transcript_id

# Column ID       Description
# ------------------------------------------------
# kgID            Known Gene ID
# mRNA            mRNA ID
# spID            UniProt protein Accession number
# spDisplayID     UniProt display ID
# geneSymbol      Gene Symbol
# refseq          RefSeq ID
# protAcc         NCBI protein Accession number
# description     Description
# rfamAcc         Rfam accession number
# tRnaName        Name from the tRNA track

sub kgXref {
    my ($xref, $key, $warn) = @_;
    $warn //= 1;
    $key //= 'mRNA';
    my @cols = qw(kgID mRNA spID spDisplayID geneSymbol refseq protAcc description rfamAcc tRnaName);
    if (!(grep { $_ eq $key } @cols)) {
        die "$key not a column name, must be one of:\n".join(@cols)
    }
    my @rest = grep { $_ ne $key } @cols;
    my $csv = Text::CSV_XS->new({sep_char   => "\t"});
    $csv->column_names(@cols);
    open(my $fh, '<', $xref) or die $!;
    my %munged;
    while (my $d = $csv->getline_hr($fh)) {
        if (! defined($d->{$key})) {
            warn "Skipping line based on empty key $key: ".join( ',', map { $d->{$_} } @rest ) if $warn;
            next;
        }
        $munged{ $d->{ $key } } = { map { $_ => $d->{$_} } @rest };
    }
    \%munged;
}

# Column ID                    Description
# ------------------------------------------------
# knownCanonical.chrom         Chromosome
# knownCanonical.chromStart    Start position (0 based). Represents transcription start for + strand genes, end for - strand genes
# knownCanonical.chromEnd      End position (non-inclusive). Represents transcription end for + strand genes, start for - strand genes
# knownCanonical.clusterId     Which cluster of transcripts this belongs to in knownIsoforms
# knownCanonical.transcript    Corresponds to knownGene name field.
# knownCanonical.protein       Accession of the associated protein, or UCSC ID in newer tables.
# kgXref.kgID                  Known Gene ID
# kgXref.mRNA                  mRNA ID
# kgXref.spID                  UniProt protein Accession number
# kgXref.spDisplayID           UniProt display ID
# kgXref.geneSymbol            Gene Symbol
# kgXref.refseq                RefSeq ID
# kgXref.protAcc               NCBI protein Accession number
# kgXref.description           Description
# kgXref.rfamAcc               Rfam accession number
# kgXref.tRnaName              Name from the tRNA track

sub knownCanonical_kgXref {
    my ($xref, $key, $warn) = @_;
    $warn //= 1;
    $key //= 'kgXref.refseq';
    my @cols = qw(knownCanonical.chrom knownCanonical.chromStart
    knownCanonical.chromEnd knownCanonical.clusterId knownCanonical.transcript
    knownCanonical.protein kgXref.kgID kgXref.mRNA kgXref.spID
    kgXref.spDisplayID kgXref.geneSymbol kgXref.refseq kgXref.protAcc
    kgXref.description kgXref.rfamAcc kgXref.tRnaName);
    if (!(grep { $_ eq $key } @cols)) {
        die "$key not a column name, must be one of:\n".join(@cols)
    }
    my @rest = grep { $_ ne $key } @cols;
    my $csv = Text::CSV_XS->new({sep_char   => "\t"});
    $csv->column_names(@cols);
    open(my $fh, '<', $xref) or die $!;
    my %munged;
    while (my $d = $csv->getline_hr($fh)) {
        if (! defined($d->{$key})) {
            warn "Skipping line based on empty key $key: ".join( ',', map { $d->{$_} } @rest ) if $warn;
            next;
        }
        $munged{ $d->{ $key } } = { map { $_ => $d->{$_} } @rest };
    }
    \%munged;
}

1;
