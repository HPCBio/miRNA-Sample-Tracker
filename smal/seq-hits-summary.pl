#!/usr/bin/env perl
use 5.014;
use strict;
use warnings;
use autodie;
use Getopt::Long;
use Text::CSV_XS;
use lib $ENV{GIT_REPO}."/hpcbio/smalheiser_scripts/v2/lib";
use miRNA_DB;
use Carp;
use List::Util qw(sum);

my $db;
my $cutoff = 0;
my $ref = 'hg19_spiked';
my $kgXref = '/home/mirrors/igenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/kgXref.txt';
my $max = 0;
my $min = 1;
my $sum = 0;

my $USAGE = <<USAGE;

    perl $0 -b <BAM> -f <ANNOTATION.bed/gff> -t <bed/gff> -o my_table.txt -d <DATABASE> -c <CUTOFF> -r <REFERENCE>

Requires -b, -f, -d, -k

Also requires modules 'perl/5.16.1_hpcbio', 'samtools', and 'bedtools2'

Works best with a queryname-sorted BAM file, which will sort the srna IDs in order

USAGE

GetOptions(
    'sum_co=s'      => \$sum,      # total count filter (collides with cutoff)
    'database=s'    => \$db,
    'cutoff=i'      => \$cutoff,  # min count for srna in any one sample
    'reference=s'   => \$ref,
    'kgXref=s'      => \$kgXref,
    'min=i'         => \$min,
    'max=i'         => \$max      # max # hits in reference; may add max # alns reported as well
);

die $USAGE if !defined($db);
die $USAGE if $sum && $cutoff; # can't define both

my $dbh = miRNA_DB::get_dbh($db);

my $nms = miRNA_DB::get_db_names($dbh);

my $buffer;

# callback to process data into chunks for this script
my $coderef = sub {
    my $st = shift;
    # get statement data
    while (1) {
        my $row = $st->fetchrow_hashref();
        # if !data,
        if (!defined $row) {
            # at end of stream so return tmp pointer to current buffer (which will be cleared out for next iteration)
            my $tmp = $buffer;
            $buffer = undef;
            return $tmp;
        } else { # else:
            # check whether there is a buffer, if so 
            # check whether buffer ID doesn't match current ID
            if ($buffer) {
                my $tmp;
                if ($buffer->{srna_id} ne $row->{srna_id}) {
                    $tmp = $buffer; # assign to old buffer
                    
                    # reassign buffer to point to new data; old data should be safe in $tmp
                    $buffer = {
                               srna_id  => $row->{srna_id},
                               count    => $row->{count},
                               sequence    => $row->{sequence}
                               };
                }                    
                $buffer->{data}{ $row->{name} } = $row->{hits} if $row->{name};
                return $tmp if $tmp;
            } else {
                $buffer = {
                           srna_id  => $row->{srna_id},
                           count    => $row->{count},
                           sequence    => $row->{sequence}
                           };
                $buffer->{data}{$row->{name}}  = $row->{hits} if $row->{name};
            }
        }
    }
};

my $it = miRNA_DB::get_db_iterator($dbh, $sum, $coderef);

say join("\t", qw(ID COUNT SEQUENCE LENGTH), @$nms, 'SUM');

my $ct = 0;
while (my $d = $it->()) {
    my @totals = map {
                exists $d->{data}{$_} ? $d->{data}{$_} : 0
                } @$nms;
    say join("\t",
             'srna'.$d->{srna_id},
             $d->{count},
             $d->{sequence},
             length($d->{sequence}),
             @totals,
             sum(@totals)
             );
}

