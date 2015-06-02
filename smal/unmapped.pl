#!/usr/bin/env perl
use 5.014;
use strict;
use warnings;
use autodie;
use Getopt::Long;
use Text::CSV_XS;
use FindBin;
use lib "$FindBin::RealBin/../lib";
use miRNA_DB;
use BamUtils;
use SampleUtils;
use Carp;
use Data::Dumper;
use List::Util qw(sum);

my $db;
my $cutoff = 0;
my $max = 0;
my $min = 1;
my $sum = 0;
my $hg19 = 'hg19_spiked';

my ($hairpin, $mature, $ac);

my ($addcts, $meta_file, $sample_order_file);

my $USAGE = <<USAGE;

Writing up usage, but here's a basic command line:

perl smalheiser_scripts/v2/smal/unmapped.pl -sum 20 \
    -d nonn.sqlite3 \
    -ac abundant_contaminant_ALL/all.name.bam \
    -add \
    -meta sample_info.txt \
    -hg19 hg19_spike

The '-hg19' option should be set to the name of the database to screen against
(e.g. we're looking for seqs unmapped to hg19, but they may be mapped to others)

Also requires modules 'perl/5.16.1_hpcbio', 'samtools', and 'bedtools2'

Works best with a queryname-sorted BAM file, which will sort the srna IDs in order

USAGE

GetOptions(
    'sum_co=s'      => \$sum,      # total count filter (collides with cutoff)
    'database=s'    => \$db,
    #'cutoff=i'      => \$cutoff,  # min count for srna in any one sample (NYI)
    'min=i'         => \$min,
    'max=i'         => \$max,      # max # hits in reference; may add max # alns reported as well
    'ac=s'          => \$ac,      # abundnant-contaminant query name-sorted BAM file
    #'mature=s'      => \$mature,   # miRBase mature query name-sorted BAM file (NYI)
    #'hairpin=s'     => \$hairpin  # miRBase hairpin query name-sorted BAM file (NYI)
    'add_cts'       => \$addcts,
    'meta=s'        => \$meta_file,
    'list=s'        => \$sample_order_file,
    'hg19=s'        => \$hg19 
);

die $USAGE if !defined($db);
die $USAGE if $sum && $cutoff; # can't define both

my $dbh = miRNA_DB::get_dbh($db);

my $nms = miRNA_DB::get_db_names($dbh);

my $buffer;

my ($ac_map, $mature_map, $hairpin_map);

if ($ac) {
    my $ids = get_seq_ids_by_total_ct($dbh, $sum);
    $ac_map = BamUtils::get_ref_map($ac, $ids);
}

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

open(my $in_mir_fh, '>', 'not-in-hg19.txt') or die $!;

my @cols = (qw(ID COUNT SEQUENCE LENGTH),
            (grep {
        $_ ne $hg19
      } @$nms),
            'abundant-contaminant hits');

if ($meta) {
    say $in_mir_fh SampleUtils::meta_data_str($meta, $sl, $#cols);
}

say $in_mir_fh join("\t", @cols, @$sl);

my $ct = 0;
while (my $d = $it->()) {
    # if it doesn't map to hg19
    if (! exists($d->{data}{$hg19}) ) {
        my $str = join("\t",
             'srna'.$d->{srna_id},
             $d->{count},
             $d->{sequence},
             length($d->{sequence}),
             (map {
                $d->{data}{$_} ? 1 : 0
                }
              grep { $_ ne $hg19 }
             @$nms),
             exists $ac_map->{ $d->{srna_id} } ? join(';', sort keys %{ $ac_map->{ $d->{srna_id} } }) : ''
             );
        if ($addcts) {
            my $counts = miRNA_DB::get_seq_counts($dbh, $d->{srna_id});
            $str .= "\t".join("\t", map { exists $counts->{$_} ? $counts->{$_} : 0 } @$sl);
        }
        say $in_mir_fh $str;
    }
}

exit;

# TODO: move to library, add cutoff/min/max filtering
sub get_seq_ids_by_total_ct {
    my ($dbh, $tc) = @_;
        
    my $mapped_cts = $dbh->prepare(<<SQL);
SELECT s.srna_id
    FROM srna AS s
    WHERE
        s.count >= ?
SQL
    my ($id);
    
    $mapped_cts->bind_columns(\$id);

    $mapped_cts->execute($tc) or die $mapped_cts->errstr;
    
    my @ids;
    while ($mapped_cts->fetch()) {
        push @ids, $id
    }
    \@ids;

}
