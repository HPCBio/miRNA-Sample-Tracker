#!/usr/bin/env perl
use 5.010;
use strict;
use warnings;
use autodie;
use File::Spec;
use File::Basename;
use Getopt::Long;
use DBI;
use Data::Dumper;
use Bio::DB::Sam;

my $USAGE = <<USAGE;

    $0 -d <DB_NAME> -c <CUTOFF>

    This script will generate a summary table of the total number of mapped hits
    per sample, taking into account the sample sequence redundancy. User must
    provide the SQLite3 database. An optional cutoff level can be added
    (default=1, or keep everything)

USAGE

{

    my $dbfile;
    my $cutoff = 1;
    my $reference = 'hg19';
    my $bam;
    my $fasta;

    GetOptions('db=s'       => \$dbfile,
               'cutoff=i'   => \$cutoff,
               'ref=s'      => \$reference,
               'bam=s'      => \$bam,
               'fasta=s'    => \$fasta);

    if (!-e $dbfile) {
        die "Database $dbfile doesn't exist";
    }

    my $seq_data = bam_info($bam, $fasta);
    
    my $dbh = DBI->connect("dbi:SQLite:$dbfile","","",{RaiseError =>1})
        or die "Couldn't connect to database: " . DBI->errstr;
    $dbh->do("PRAGMA foreign_keys = ON");
    
    $dbh->do('PRAGMA synchronous = 0');      # Non transaction safe!!!
    $dbh->do('PRAGMA cache_size = 8000000'); # 2 GB dynamic cache increase
                                             # makes index creation faster
    
    mapped_seqs($dbh, $cutoff, $reference, $seq_data);
}

sub mapped_seqs {
    my ($dbh, $cutoff, $reference, $seq_data) = @_;

    # get sequences from specified database based on minimal count
    
    my $dbid = get_db_names($dbh, $reference);
    
    my $mapped_cts = $dbh->prepare(<<SQL);
SELECT DISTINCT s.srna_id, s.sequence
    FROM srna AS s
    JOIN sample2srna AS sm2s ON (s.srna_id=sm2s.srna_id)
    JOIN srna2db AS s2d ON (s.srna_id=s2d.srna_id)
    WHERE
        s2d.db_id = ?
    AND
        sm2s.count >= ?
    ORDER BY s.srna_id
SQL

    my ($id, $count, $seq);
    
    $mapped_cts->bind_columns(\$id, \$seq);

    $mapped_cts->execute($dbid, $cutoff) or die $mapped_cts->errstr;
    
    while ($mapped_cts->fetch()) {
        say ">srna$id ".join(',', @{$seq_data->{"srna$id"}});
        say $seq;
        #say join("\t", $id, $count, $seq);
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

sub bam_info {
    my ($bam, $fasta) = @_;
    my $db = Bio::DB::Bam->open($bam);
    
    my $h = $db->header;
    my $targets = $h->target_name;
    my $index = Bio::DB::Bam->index_open($bam);
    
    my %seq_data;
    
    my $cb = sub {
        my $alignment = shift;
        my $nm = $alignment->qname;
        my $seqid       = $targets->[$alignment->tid];
        push @{$seq_data{$nm}}, $seqid;
        #print $alignment->qname," aligns to $seqid\n";
    };
    
    for my $id (@$targets) {
        $index->fetch($db, $h->parse_region($id), $cb);
        say STDERR "Finished $id";
    }
    
    \%seq_data;
}
