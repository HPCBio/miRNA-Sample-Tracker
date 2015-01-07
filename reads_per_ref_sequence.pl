#!/usr/bin/env perl
use 5.010;
use strict;
use warnings;
use autodie;
use File::Spec;
use File::Basename;
use Getopt::Long;
use DBI;
use Bio::DB::Sam;
use Data::Dumper;
use List::MoreUtils qw(uniq);

my $USAGE = <<USAGE;

    $0 -d <DB_NAME> -r <REFERENCE>

    This script will generate a summary table of the total number of mapped hits
    per sample, taking into account the sample sequence redundancy. User must
    provide a BAM file and the SQLite3 database.
    
    NOTE: reads could technically align to multiple reference sequences. This
    does note take into account multi-mappiness, just whether a sequence has
    aligned to a reference sequence and how this translates to counts in the
    original sample
    
    
    Rows are sample IDs, columns are reference sequence IDs.

USAGE

{

    my ($dbfile, $bam);
    
    # this is the lower count limit, keep anything where an sRNA sequence is present
    # greater than or equal to this in any one sample
    my $cutoff = 1;

    GetOptions('db=s'       => \$dbfile,
               'cutoff=i'    => \$cutoff,
               'bam=s'      => \$bam);

    if (!-e $dbfile) {
        die "Database $dbfile doesn't exist";
    }

    if (!-e $bam) {
        die "Must provide BAM file"
    }
    
    # init SQLite db
    my $dbh = DBI->connect("dbi:SQLite:$dbfile","","",{RaiseError =>1})
        or die "Couldn't connect to database: " . DBI->errstr;
    $dbh->do("PRAGMA foreign_keys = ON");

    $dbh->do('PRAGMA synchronous = 0');      # Non transaction safe!!!
    $dbh->do('PRAGMA cache_size = 8000000'); # 2 GB dynamic cache increase
    
    # init BAM db
    my $bam_b = Bio::DB::Bam->open($bam);
    my $bam_i = Bio::DB::Bam->index_open($bam);
    
    my $bam_h = $bam_b->header;
    
    # iterate through ref sequences and grab unique sequences IDs for each
    my @refids = @{$bam_h->target_name};

    # retrieve counts for each; note the callback
    my $cb = sub {
        my ($align,$data) = @_;
        my $ids = $data->[0];
        push @{$ids}, $align->qname;
    };
    
    # hold the summary data
    my %summary;
    
    for my $tid (0..$#refids) {
        my $seqid = $refids[$tid];
        my $len = $bam_h->target_len->[$tid];
        
        my @ids;
        
        # grab reads from each region
        my $code = $bam_i->fetch($bam_b, $tid, 0, $len, $cb, [\@ids]);
    
        #say Dumper [uniq @ids];
        # grab counts for all reads along with
        
        # TODO: getting ugly...
        get_srna2sample($dbh, [ uniq @ids ], $seqid, \%summary, $cutoff);
    }
    
    open(my $outfh, '>', "$bam.sample_summary.$cutoff.txt") or die $!;
    
    say $outfh join("\t", "Sample", @refids);
    
    for my $smid (sort keys %summary) {
        say $outfh join("\t", $smid, map { exists $summary{$smid}{$_} ? $summary{$smid}{$_} : 0 } @refids);
    }
}

sub get_srna2sample {
    my ($dbh, $ids, $seqid, $summary, $cutoff) = @_;
    
    my $list = join(',', (map { s/srna//; $_ } @$ids));
    
    my $sql = <<SQL;
SELECT sm.name, sm.sample_id, SUM(sm2s.count) AS Total
    FROM sample AS sm
    JOIN sample2srna AS sm2s ON (sm.sample_id=sm2s.sample_id)
    WHERE
        sm2s.count >= ?
    AND
        sm2s.srna_id IN (
            $list
        )
    GROUP BY sm.sample_id
SQL
    
    my $mapped_cts = $dbh->prepare($sql);
    
    my ($nm, $smid, $total);
    
    $mapped_cts->bind_columns(\$nm, \$smid, \$total);
    
    $mapped_cts->execute($cutoff);
    
    while ($mapped_cts->fetch()) {
        $summary->{$nm}->{$seqid} = $total;
    }
    
    say STDERR "Finished $seqid";
    
}
