#!/usr/bin/env perl

use 5.014;
use strict;
use warnings;
use autodie;
use Data::Dumper;
use Getopt::Long;
use PerlIO::gzip;
no if $] >= 5.017011, warnings => 'experimental::smartmatch';
use lib $ENV{GIT_REPO}."/hpcbio/smalheiser_scripts/v2/lib";
use miRNA_DB;
use BamUtils;

my $max_alns = my $min_alns = my $max_rpts = 0;
my $ref = 'unknown';
my $cutoff;
my $dbfile;
my $bam;
my $sum = 0;

GetOptions(
    'bam=s'             => \$bam,
    'db=s'              => \$dbfile,
    'ref=s'             => \$ref,
    'cutoff=i'          => \$cutoff,
    'sum=i'             => \$sum,
    'max_alns=i'        => \$max_alns,
    'min_alns=i'        => \$max_alns,
    'max_reports=i'     => \$max_rpts
);

my $USAGE = <<USAGE;

    $0 -d <DB_NAME> -c <CUTOFF> -b <BAM> -r <REF>

    This script will tak an arbitrary BAM file and output a table of the hits to
    the reference, including basic information such as the number of hits in
    that reference, location/strand, sequence, type of alignment, etc.
    
    Note this requires a database reference as well (-r).  This requirement will
    be removed at a future point, as it is essentially redundant with the BAM file
    name.

USAGE

my $dbh = miRNA_DB::get_dbh($dbfile);

# grab IDs with minimal cutoff from database
my %ids;

if ($sum || $cutoff) {
    %ids = $sum ? map { $_ => 1 } @{ miRNA_DB::get_seq_ids_by_total_count($dbh, $ref, $sum) } :
                  map { $_ => 1 } @{ miRNA_DB::get_seq_ids($dbh, $ref, $cutoff) } ;
}

# open BAM stream and skip unmapped reads
open(my $bamfh, "samtools view -F 4 $bam |");

open(my $logfh, "| gzip > $ref.tbl.gz");

say $logfh join("\t", qw|READ DB SEQUENCE LENGTH STRAND MATCH_TYPE NUM_ALNS REFERENCE START END MATCH_STRING(CIGAR)
     HAIRPIN_POS HAIRPIN_SCORE |);

my %rpt_cache;

SAM_STREAM:
while (my $ln = <$bamfh>) {
    if ($ln =~ /^@/) {
        next;
    }
    my @data = split("\t",$ln,12);
    
    (my $id = $data[0]) =~ s/^srna//;
    
    if ($sum || $cutoff) {
        next unless exists $ids{$id};
    }
    
    # note this is very much for SE data only :P
    my $strand = BamUtils::flag_strand($data[1]);
    my $type = BamUtils::cigar_stats($data[5]);
    my $len = BamUtils::cigar_len($data[5]);
    my $tags = BamUtils::tag_info($data[11]);
    
    my $end =  $data[3] + $len - 1;
    if ($max_alns && ($tags->{alns} > $max_alns)) {
        next;
    }

    #if ($min_alns && $tags->{alns} < $min_alns) {
    #    next;
    #}
    
    if ($max_rpts && ($tags->{alns} >= $max_rpts)) {
        next if $rpt_cache{$data[0]}++ >= $max_rpts;
    }
    
    say $logfh join("\t",
        $data[0],
        $ref,
        $data[9] || '',
        length($data[9]),
        $strand,
        $type,
        $tags->{alns},
        @data[2, 3],
        $end,
        $data[5],
        $tags->{hploc} // 0,
        $tags->{hpscore} // 0
        );
}

