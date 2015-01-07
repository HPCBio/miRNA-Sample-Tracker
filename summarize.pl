#!/usr/bin/env perl

use 5.014;
use strict;
use warnings;
use autodie;
use Data::Dumper;
use Getopt::Long;
use PerlIO::gzip;
use Bio::DB::Fasta;
no if $] >= 5.017011, warnings => 'experimental::smartmatch';
use lib $ENV{GIT_REPO}."/hpcbio/smalheiser_scripts/v2/lib";
use miRNA_DB;

my $max_alns = 0;
my $dbfile;
my $bam;
my $ref;
my $cutoff;

GetOptions(
    'bam=s'         => \$bam,
    'db=s'          => \$dbfile,
    'ref=s'         => \$ref,
    'cutoff=i'      => \$cutoff
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
my %ids = map { $_ => 1 } @{ miRNA_DB::get_seq_ids($dbh, $ref, $cutoff) };

open(my $bamfh, "samtools view $bam |");

open(my $logfh, "| gzip > $ref.tbl.gz");

say $logfh join("\t", qw|READ DB SEQUENCE LENGTH STRAND MATCH_TYPE NUM_ALNS REFERENCE START END MATCH_STRING(CIGAR)
     HAIRPIN_POS HAIRPIN_SCORE |);

SAM_STREAM:
while (<$bamfh>) {
    if (/^@/) {
        next;
    }
    my @data = split("\t",$_,12);
    
    (my $id = $data[0]) =~ s/^srna//;
    
    next unless exists $ids{$id};

    my ($type, $strand);

    for my $l ($data[1]) {
        if ($l == 4) {
            next SAM_STREAM
        } elsif ($l == 0) {
            $strand = '+'
        } elsif ($l == 256) {
            $strand = '+'
        } elsif ($l == 16) {
            $strand = '-'
        } elsif ($l == 272) {
            $strand = '-'
        }
    }

    die "Line $_" unless $strand;

    # exact match, check CIGAR string
    if ($data[5] =~ /^\d+M$/) {
        #print $exactfh $_;
        $type = 'EXACT';
    } elsif ($data[5] =~ /^(\d+[HS])?\d+M(\d+[HS])?$/) {
        #print $exactcfh $_;
        $type = 'CLIPPED';
    } else {
        #print $mmfh $_;
        $type = 'MISMATCH/INDEL';
    }
    
    # TODO: should pull this out into a method!!!!
    my $len = 0;
    while ($data[5] =~ /(\d+)[MDN=XP]/g) {
        $len += $1;
    }

    my ($alns, $hpscore, $hploc) = (1,0,0);

    while ($data[11] =~ /(Z[NHL]):\w:(\S+)/g) {
        for ($1) {
            when ('ZH') {
                $hpscore = $2;
            }
            when ('ZN') {
                $alns = $2;
            }
            when ('ZL') {
                $hploc = $2;
            }
        }
    }
    my $end =  $data[3] + $len - 1;
    say $logfh join("\t", $data[0], $ref, $data[9] || '', length($data[9]), $strand, $type, $alns, @data[2, 3], $end ,$data[5], $hploc, $hpscore);
}

#sub get_seq_ids {
#    my ($dbh, $reference, $cutoff) = @_;
#
#    # get sequences from specified database based on minimal count
#    
#    my $dbid = get_db_names($dbh, $reference);
#    
#    # TODO: optimize into a column call
#    my $mapped_cts = $dbh->prepare(<<SQL);
#SELECT DISTINCT s.srna_id
#    FROM srna AS s
#    JOIN sample2srna AS sm2s ON (s.srna_id=sm2s.srna_id)
#    JOIN srna2db AS s2d ON (s.srna_id=s2d.srna_id)
#    WHERE
#        s2d.db_id = ?
#    AND
#        sm2s.count >= ?
#SQL
#
#    my ($id, $count, $seq);
#    
#    $mapped_cts->bind_columns(\$id);
#
#    $mapped_cts->execute($dbid, $cutoff) or die $mapped_cts->errstr;
#    
#    my @ids;
#    while ($mapped_cts->fetch()) {
#        push @ids, $id
#    }
#    \@ids;
#}


#sub get_db_names {
#    my ($dbh, $ref) = @_;
#    my $id = $dbh->selectcol_arrayref(<<SQL);
#SELECT db_id FROM db_info
#WHERE
#    db_info.name="$ref"
#SQL
#    die "No ID returned for $ref" if !defined($id) || @$id == 0;
#    @$id[0];
#}