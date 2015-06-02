#!/usr/bin/env perl
use 5.010;
use strict;
use warnings;
use autodie;

use FindBin;
use lib "$FindBin::Bin/../lib";

use File::Spec;
use File::Basename;
use Getopt::Long;
use DBI;

my $USAGE = <<USAGE;

    $0 -f <sample.compressed.fasta> -db <DB_NAME> --dir <DIR> 

This reads a FASTA from a collapsed sample sequence file and a database
(created using db-init.pl) and loads a table with the list of read ID
and counts for that sample. Match is based upon the sequence read.

If given a directory this will simply iterate through all the files and
(sequentially) determine counts for each.  It ridiculously NON-parallel at this
time.

TODO: optimize this into a table for the database?

USAGE


{

    my $dir;
    my $dbfile;
    my $ext = '.derep.fasta';

    GetOptions('db=s'       => \$dbfile,
               'dir=s'      => \$dir,
               'ext=s'      => \$ext);

    die unless (defined $dir && -e $dir);

    # TODO: allow other extensions?
    my @fastas = glob(File::Spec->catfile($dir, '*.derep.fasta'));

    die "Couldn't find fasta files using $dir with file extension $ext"
        unless @fastas;

    if (!-e $dbfile) {
        die "Database $dbfile doesn't exist";
    }

    my $dbh = DBI->connect("dbi:SQLite:$dbfile","","",{RaiseError =>1})
        or die "Couldn't connect to database: " . DBI->errstr;
    $dbh->do("PRAGMA foreign_keys = ON");

    $dbh->do('PRAGMA synchronous = 0');      # Non transaction safe!!!
    $dbh->do('PRAGMA cache_size = 4000000'); # 2 GB dynamic cache increase
                                             # makes index creation faster

    for my $f (@fastas) {
        say STDERR "Loading data from file $f";
        $dbh->do("BEGIN");
        load_counts($dbh, $f);
        $dbh->do("COMMIT");
    }

    $dbh->do("BEGIN");
    $dbh->do ("CREATE INDEX sample2srna_sample_idx on sample2srna( sample_id )");
    $dbh->do("COMMIT");

}

sub load_counts {
    my ($dbh, $fasta_file) = @_;
    my $bn = File::Basename::basename($fasta_file);
    $bn =~ s/_[ATCG]{6}_L\d+_.*//;

    my $sample_id = add_sample($dbh, $bn);

    my $get_seq = $dbh->prepare(<<SQL);
SELECT srna_id
    FROM srna
    WHERE sequence = ?;
SQL

    my $add_counts = $dbh->prepare(<<SQL);
INSERT INTO sample2srna (srna_id, sample_id, count)
    VALUES(?, ?, ?)
SQL
    open(my $FASTA, '<', $fasta_file);
    
    my ($srna_id);

    $get_seq->bind_columns(\$srna_id);
    while (!eof($FASTA)) {
        my $header = <$FASTA>;
        my $seq = <$FASTA>;
        my ($name, $count);
        if ($header =~ /^>(\d+)-(\d+)$/) {
            ($name, $count) = ($1, $2);
        } else {
            die "Sequence header did not match: $header";
        }
        chomp $seq;
        $get_seq->execute($seq) or die $get_seq->errstr;
        $get_seq->fetch() or die $get_seq->errstr;
        $add_counts->execute($srna_id, $sample_id, $count) or die $add_counts->errstr;
    }
    close $FASTA;
    1;
}

sub add_sample {
    my ($dbh, $nm) = @_;
    my $add_sample = $dbh->prepare(<<SQL);
INSERT INTO sample ( name )
    VALUES(?)
SQL
    $add_sample->execute($nm);
    my $id = $dbh->last_insert_id('','','','');
    $id;
}
