#!/usr/bin/env perl
use 5.010;
use strict;
use warnings;
use autodie;
use File::Spec;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(sum);
use lib $ENV{GIT_REPO}."/hpcbio/smalheiser_scripts/v2/lib";

use FindBin;
use lib "$FindBin::Bin/../lib";

use miRNA_DB;

my $USAGE = <<USAGE;

    $0 -d <DB_NAME> -c <CUTOFF> -r <REFERENCE> -u 

This script will generate a count table (# times a particular small RNA is
present in a specific sample).

The following parameters are accepted:
    
    -c|--cutoff   : a minimal cutoff, which represents the minimal count found in any one sample
    --sum         : sum cutoff, total counts for sequence in all samples
    -r|--ref      : reference database to limit sequence data from
    -u|--unmapped : unmapped only (reference database is ignored if this is set)
    --sample      : tab-delimited file with sample information
    -b|--by       : name of column to sort by in sample table
    -l|--list     : list of sample names to sort by (overrides -s and -b)

USAGE

{

    my $dbfile;
    my $list;
    #my $outfile = 'counts.txt';
    my $cutoff = 1;
    my $ref;
    my $unmapped = 0;
    my $sample;
    my $by;

    GetOptions('db=s'       => \$dbfile,
               'cutoff=i'   => \$cutoff,
               'ref=s'      => \$ref,
               'unmapped'   => \$unmapped,
               'sample=s'   => \$sample,
               'by=s'       => \$by,
               'list=s'     => \$list);

    my $dbh = miRNA_DB::get_dbh($dbfile);
    
    my ($sample_data, $sample_cgs, $sample_names);
    
    if ($list) {
        # custom list, retaining order
        $sample_names = miRNA_DB::get_sample_list($list);
    } else {
        # grab from the database
        $sample_names = miRNA_DB::get_sample_names($dbh);        
    }
    
    if ($sample) {
        ($sample_data, $sample_cgs) = miRNA_DB::get_sample_info($sample);
        if ($by && (! grep { $_ eq $by } @$sample_cgs )) {
            die "Provided sort category (via -b/--by) not in sample table"
        }
    }
        
    my $dbid;
    
    if ($ref) {
        $dbid = miRNA_DB::get_dbid($dbh, $ref);
    }
    
    my $get_srnaids;
    
    my ($srnaid, $srnaid_seq, $srnaid_ct, $sample_num, $db_hits);
    
    if ($unmapped) {
        
        $get_srnaids = $dbh->prepare(<<SQL1);

SELECT s.srna_id, s.sequence, s.count, COUNT(*) as sample_num
FROM srna as s
JOIN sample2srna as s2r ON (s.srna_id=s2r.srna_id)
WHERE
    NOT EXISTS (
        SELECT * FROM srna2db AS r2d WHERE r2d.srna_id=s.srna_id
    )
    AND
    s2r.count >= ?
GROUP BY s.srna_id
ORDER BY s.srna_id

SQL1

        $get_srnaids->bind_columns(\$srnaid, \$srnaid_seq, \$srnaid_ct, \$sample_num);

    } elsif (defined($dbid)) {
        
        $get_srnaids = $dbh->prepare(<<SQL2);
        
SELECT s.srna_id, s.sequence, s.count, COUNT(*) as sample_num, r2d.hits
FROM srna as s
JOIN sample2srna as s2r ON (s.srna_id=s2r.srna_id)
JOIN srna2db as r2d ON (s.srna_id=r2d.srna_id)
WHERE
    r2d.db_id = $dbid
    AND
    s2r.count >= ?
GROUP BY s.srna_id
ORDER BY s.srna_id
SQL2

        $get_srnaids->bind_columns(\$srnaid, \$srnaid_seq, \$srnaid_ct, \$sample_num, \$db_hits);

    } else {
        die "Undefined database ID"
    }

    my $get_cts = $dbh->prepare(<<SQL);
SELECT sm.name, sm2r.count
    FROM sample2srna sm2r
    JOIN sample sm ON (sm2r.sample_id=sm.sample_id)
    WHERE sm2r.srna_id=?
SQL
    
    $get_srnaids->execute($cutoff) or die $get_srnaids->errstr;

    if (defined $sample_data) {
        if ( defined $by ) {
            if ($by eq 'Sample_ID' ) {
                # sort by ID
                @{$sample_names} = sort {
                    $a <=> $b
                } @{$sample_names};
            } else {
                # sort by custom column
                @{$sample_names} = sort {
                    $sample_data->{$a}{$by} cmp $sample_data->{$b}{$by}
                    ||
                    $a <=> $b
                } @{$sample_names};
                say join("\t", '', '', '', '', $by, map { $sample_data->{$_}{$by} || '' } @{ $sample_names }, '');
            }
        } else {
            # this outputs all the sample meta data above the columns and allows for a custom column sort
            for my $cat (@$sample_cgs) {
                if ($cat eq 'Sample_ID') {
                    say join("\t", '', '', '', '', $cat, @{ $sample_names }, '');
                } else {
                    say join("\t", '', '', '', '', $cat, map { $sample_data->{$_}{$cat} || '' } @{ $sample_names }, '');
                }
            }
        }
    }
    
    say join("\t", "small_RNA_ID", 'small_RNA_Seq', 'Length', 'Num_Samples', 'DB Hits', @{ $sample_names }, 'Sum_Total');
    
    while ($get_srnaids->fetch()) {
        $get_cts->execute($srnaid) or die $get_cts->errstr;
        my $d = $get_cts->fetchall_hashref('name');
        my @vals = map { exists $d->{$_} ? $d->{$_}->{count} : 0 } @{ $sample_names };
        say join("\t",
                 "srna$srnaid",
                 $srnaid_seq,
                 length($srnaid_seq),
                 $sample_num,
                 $db_hits // 0,
                 @vals,
                 (sum @vals));
    }
}
