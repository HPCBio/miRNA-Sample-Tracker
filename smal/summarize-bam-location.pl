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
use SampleUtils;

my $max_alns = my $min_alns = my $max_rpts = 0;
my $ref = 'unknown';
my $cutoff;
my $dbfile;
my $bam;
my $sum = 0;
my $addcts = 0;
my $meta_file;
my $sample_order_file;

GetOptions(
    'bam=s'             => \$bam,
    'db=s'              => \$dbfile,
    'ref=s'             => \$ref,
    'cutoff=i'          => \$cutoff,
    'sum=i'             => \$sum,
    'max_alns=i'        => \$max_alns,
    'min_alns=i'        => \$max_alns,
    'max_reports=i'     => \$max_rpts,
    'add_cts'           => \$addcts,
    'meta=s'            => \$meta_file,
    'list=s'            => \$sample_order_file
);

my $USAGE = <<USAGE;

    $0 -d <DB_NAME> -c <CUTOFF> -b <BAM> -r <REF>

This script will take an arbitrary query-sorted BAM file and output a table of
the hits to the reference, including basic information such as the ID and
sequence, number of hits in the reference, and location/strand.  At this time
there is no check performed on the input file, but logic in this script requires
sorting in this order.

Only one line is reported for each small RNA ID. Therefore it's highly
recommended to filter the data under instances where there are expected to be
multiple alignments to reference databases. This can be done using:

* 'cutoff' - minimal # sequences in any one sample (default = 0, or report
  everything)
* 'sum' - overall # sequences in all samples (default = 0, or report everything)
* 'max_alns' - report sequences with total # alignments less than or equal to
  this number (default = 0, or report everything)
* 'min_alns' - report sequences with total # alignments greater than or equal to
  this number (default = 0, or report everything)
* 'max_reports' (NYI) - only report this many alignments (default = 0, or report
  everything)

'cutoff' and 'sum' cannot be used together for obvious reasons, but all others
can be combined.  'max_reports' is not yet implemented but will be added.

'add_cts' can be used to add read count information to the file as well.

Note this requires a database reference as well (-r). This requirement will be
removed at a future point, as it is essentially redundant with the BAM file
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
my %rpt_cache;

my $ct = 0;

my $buffer;

my $sl = [];

my ($meta, $metacols);

if ($addcts) {
    if ($sample_order_file) {
        $sl = miRNA_DB::get_sample_list($dbh, $sample_order_file);
    } else {
        $sl = miRNA_DB::get_sample_names($dbh);
    }
    
    if ($meta_file) {
        ($meta, $metacols) = SampleUtils::get_sample_info($meta_file);
    }
}

open(my $bamfh, "samtools view -F 4 $bam |");

open(my $logfh, "| gzip > $ref.tbl.gz");

if ($meta) {
    say $logfh SampleUtils::meta_data_str($meta, $sl, 4);
}

say $logfh join("\t", qw|READ SEQUENCE LENGTH NUM_ALNS LOCATION|, @$sl);

my $coderef = sub {
    my $st = shift;
    # get statement data
    while (1) {
        my $row = $st->();
                
        # if !data,
        if (!defined $row) {
            # at end of stream so return tmp pointer to current buffer (which will be cleared out for next iteration)
            my $tmp = $buffer;
            $buffer = undef;
            return $tmp;
        } else { # else:
            # check whether there is a buffer, if so 
            # check whether buffer ID doesn't match current ID
            next if $row->[1] & 0x4; # unmapped
            if ($buffer) {
                my $tmp;
                if ($buffer->{srna_id} ne $row->[0]) {
                    $tmp = $buffer; # assign to old buffer
                    
                    # reassign buffer to point to new data; old data should be safe in $tmp
                    $buffer = {
                               srna_id  => $row->[0],
                               seq      => $row->[9],
                               hits     => exists($row->[11]->{NH}) ? $row->[11]->{NH}{data} : 1
                               };
                }                    
                push @{$buffer->{data}}, $row;
                return $tmp if $tmp;
            } else {
                $buffer = {
                           srna_id  => $row->[0],
                           seq      => $row->[9],
                           hits     => exists($row->[11]->{NH}) ? $row->[11]->{NH}{data} : 1
                           };
                push @{$buffer->{data}}, $row;
            }
        }
    }
};

# TODO: could be much faster than this implementation if using Bio-Samtools, but this pure-perl version works
my $it = BamUtils::get_bam_iterator($bam, $coderef);

while (my $d = $it->()) {
    
    (my $id = $d->{srna_id}) =~ s/^srna//;
    
    if ($sum || $cutoff) {
        next unless exists $ids{$id};
    }
    if ($max_alns && ($d->{hits} > $max_alns)) {
        delete $ids{ $id };
        next;
    }
    
    if ($min_alns && $d->{hits} < $min_alns) {
        delete $ids{ $id };
        next;
    }
    #if ($max_rpts && ($tags->{hits} >= $max_rpts)) {
    #    next if $rpt_cache{$data[0]}++ >= $max_rpts;
    #}

    my @locs;
    for my $hit (@{$d->{data}}) {
        # $hit->[5] is CIGAR str
        my $len = BamUtils::cigar_len($hit->[5]);
        # length = end - start + 1, so ...
        my $end = $hit->[3] + $len - 1;
        # if this flag is true, then the seq is rev-comped for the match
        my $strand = ($hit->[1] & 0x10) ? -1 : 1;
        push @locs, sprintf("%s:%d-%d(%d)", @{$hit}[2, 3], $end, $strand)
    }
    
    my $seq_data = miRNA_DB::get_seq_info($dbh, $id, $ref);
    
    #say Dumper $seq_data;
        
    my $str = join("\t", $d->{srna_id},
                   $seq_data->{sequence},
                   length($seq_data->{sequence}),
                   $d->{hits},
                   join(';', @locs)
                   );
    if ($addcts) {
        my $counts = miRNA_DB::get_seq_counts($dbh, $id, $ref);
        $str .= "\t".join("\t", map { exists $counts->{$_} ? $counts->{$_} : 0 } @$sl);
    }
    say $logfh $str;
    say STDERR scalar(keys %ids);
    
    delete $ids{ $id };
    last if keys %ids == 0;
}
