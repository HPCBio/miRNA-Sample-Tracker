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

    $0 -d <DB_NAME> -r <REFERENCE> -c <CUTOFF>

    This script will print a table of sequences along with a list of targets IDs
    mapped to. It requires a BAM file and the SQLite database.  An optional
    cutoff can also be provided for filtering low-count sequences.

    TODO: maybe include a hit cutoff as well

USAGE

{

    my ($dbfile, $bam, $ref);
    
    my $maxhits = 0;
    
    # this is the lower count limit, keep anything where an sRNA sequence is present
    # greater than or equal to this in any one sample
    my $cutoff = 1;

    GetOptions('db=s'       => \$dbfile,
               'ref=s'      => \$ref,
               'cutoff=i'    => \$cutoff,
               'maxhits=i'     => \$maxhits,
               'bam=s'      => \$bam);

    if (!-e $dbfile) {
        die "Database $dbfile doesn't exist";
    }

    if (! defined $ref) {
        # TODO: this is a bit of a hack; preferentially we would have the file
        # name stored in the db_info table but it isn't yet, so we call this
        # explicitly for now
        
        die "Must provide reference name";
    }
    
    if (!-e $bam) {
        die "Must provide BAM file"
    }
    
    # init SQLite db
    my $dbh = DBI->connect("dbi:SQLite:$dbfile","","",{RaiseError =>1})
        or die "Couldn't connect to database: " . DBI->errstr;
    $dbh->do("PRAGMA foreign_keys = ON");

    $dbh->do('PRAGMA synchronous = 0');      # Non transaction safe!!!
    $dbh->do('PRAGMA cache_size = 8000000'); # 8 GB dynamic cache increase

    my %summary;

    get_srna_seq($dbh, \%summary, $ref, $cutoff, $maxhits);
    
    # init BAM db
    my $bam_b = Bio::DB::Bam->open($bam);
    my $bam_i = Bio::DB::Bam->index_open($bam);
    
    my $bam_h = $bam_b->header;
    
    # iterate through ref sequences and grab unique sequences IDs for each
    my @refids = @{$bam_h->target_name};

    # hold the summary data
    open(my $outfh, '>', "$bam.seqtable.$cutoff.txt") or die $!;
    say $outfh join("\t", "SRNA ID", "Sequence", "Reference ID", "Strand", "Num Hits");
    
    # retrieve counts for each; note the callback
    my $cb = sub {
        my ($align,$data) = @_;
        my $sum = $data->[0];
        my $qn = $align->qname;
        if (exists $sum->{$qn}) {
            my $f = $align->flag;
            my $strand = $f & 16 ? '-' : '+';
            say $outfh join("\t",
                            $qn,
                            $summary{$qn}{seq},
                            $data->[1][$align->tid],
                            $strand,
                            $summary{$qn}{hits});
        }
    };
    
    for my $tid (0..$#refids) {
        my $seqid = $refids[$tid];
        my $len = $bam_h->target_len->[$tid];
        
        # grab reads from each region
        my $code = $bam_i->fetch($bam_b, $tid, 0, $len, $cb, [ \%summary, \@refids ]);
    }
}

sub get_srna_seq {
    my ($dbh, $summary, $reference, $cutoff, $maxhits) = @_;
    
    my $dbid = get_db_names($dbh, $reference);

    my $sql;
    
    if ($maxhits == 0) {    
    
        $sql = <<SQL;
SELECT s.srna_id, s.sequence, s2d.hits
    FROM srna AS s
    JOIN sample2srna AS sm2s ON (s.srna_id=sm2s.srna_id)
    JOIN srna2db AS s2d ON (s.srna_id=s2d.srna_id)
    WHERE
        s2d.db_id = ?
    AND
        sm2s.count >= ?
SQL
    } else {
        
        $sql = <<SQL;
SELECT s.srna_id, s.sequence, s2d.hits
    FROM srna AS s
    JOIN sample2srna AS sm2s ON (s.srna_id=sm2s.srna_id)
    JOIN srna2db AS s2d ON (s.srna_id=s2d.srna_id)
    WHERE
        s2d.db_id = ?
    AND
        s2d.hits <= ?
    AND
        sm2s.count >= ?
SQL
    }
    
    my $mapped_cts = $dbh->prepare($sql);
    
    my ($sid, $seq, $hits);
    
    $mapped_cts->bind_columns(\$sid, \$seq, \$hits);
    
    my @args = ($dbid, $cutoff);
    
    if ($maxhits != 0) {
        unshift @args, $maxhits;
    }
    
    $mapped_cts->execute(@args);
    
    while ($mapped_cts->fetch()) {
        $summary->{"srna$sid"}{seq} = $seq;
        $summary->{"srna$sid"}{hits} = $hits;
    }
}

sub get_db_names {
    my ($dbh, $ref) = @_;
    my $id = $dbh->selectcol_arrayref(<<SQL);
SELECT db_id FROM db_info
WHERE
    db_info.name="$ref"
SQL
    die "No ID returned for $ref" if !defined($id) || @$id == 0;
    @$id[0];
}

