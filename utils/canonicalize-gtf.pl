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

my $USAGE = <<USAGE;

    perl $0 -f <ANNOTATION.bed/gff> -k <KNOWNCANONICAL TABLE>

Requires -b, -k

Reads a GTF file and writes (to STDOUT) a canonicalized GTF for use in
downstream analyses.  Requires a specialized mapping of the knownCanonical
data to the kgXref tables for UCSC hg19.

USAGE

GetOptions(
    'features=s'    => \$features,
    'knownCanonical=s' => \$kc,
    'exon'          => \$exon
);

open( my $fh, '<', $features) or die $!;

my $canonical = FeatureUtils::knownCanonical_kgXref($kc);

my $ct = 0;
while (<$fh>) {
    chomp;
    my @data = split("\t", $_);
    next if $exon && $data[2] ne 'exon';
    my $atts = FeatureUtils::gtf_attributes($data[8]);
    if (!exists $atts->{transcript_id}) {
        die "Missing transcript_id for $_";
    }
    if (exists($canonical->{ $atts->{transcript_id }[0]} )) {
        say join("\t", @data);
    }
}

