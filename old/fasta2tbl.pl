#!/usr/bin/env perl
use 5.010;
use strict;
use warnings;
use autodie;
use File::Spec;
use Getopt::Long;
use DBI;

{

my $USAGE = <<USAGE;

    $0 -f <final.compressed.fasta> -d <DB_NAME> -c

This script is used to create a table from the collapsed FASTQ for bulk-loading
into a SQLite database; this will be used to deconvolute counts per sample.

TODO: optimize this into

USAGE

my $fasta;
my $dbfile = 'srna.sqlite3';
my $create = 0;

GetOptions('fasta=s'    => \$fasta,
           'db=s'       => \$dbfile,
           'create'     => \$create);

die unless defined $fasta && -e $fasta;

if (-e $dbfile) {
    die "Database $dbfile exists, set --create to overwrite" if !$create;
    unlink $dbfile;
}

print "Creating new DNA Database\n";

my $dbh = DBI->connect("dbi:SQLite:$dbfile","","",{RaiseError =>1})
    or die "Couldn't connect to database: " . DBI->errstr;

$dbh->do (<<SQL);
CREATE TABLE srna
    (
    id        integer(10) PRIMARY KEY AUTOINCREMENT,
    name      varchar(15),
    count     integer(10),
    sequence  varchar(50)
    );

CREATE TABLE sample
    (
    id        integer(10) PRIMARY KEY AUTOINCREMENT,
    name      varchar(15)
    );

CREATE TABLE sample2srna
    (
    srna_id      integer(10),
    sample_id    integer(10),
    count        integer(10)
    );

CREATE TABLE db_info
    (
    id        integer(10) PRIMARY KEY AUTOINCREMENT,
    name      varchar(15)
    );

CREATE TABLE srna2db
    (
    srna_id      integer(10),
    db_id        integer(10),
    hits         integer(10)
    );

SQL

$dbh->do('PRAGMA synchronous = 0');      # Non transaction safe!!!
$dbh->do('PRAGMA cache_size = 2000000'); # 2 GB dynamic cache increase
                                         # makes index creation faster

say STDERR "Loading small RNA data";
$dbh->do("BEGIN");
import_srna($dbh, $fasta);
$dbh->do("COMMIT");

say STDERR "Indexing small RNA data";
$dbh->do("BEGIN");
$dbh->do ("CREATE UNIQUE INDEX srna_seq_idx ON srna (sequence)");
$dbh->do ("CREATE UNIQUE INDEX srna_name_idx ON srna (name)");
$dbh->do("COMMIT");

}

sub import_srna {
    my ($dbh, $fasta_file) = @_;
    my $add = $dbh->prepare("INSERT INTO srna ( name, count, sequence ) VALUES(?,?,?)");

    open(my $FASTA, '<', $fasta_file);

    while (!eof($FASTA)) {
        my $header = <$FASTA>;
        my $seq = <$FASTA>;
        my ($name, $count);
        if ($header =~ /^>(\S+)\s(\d+)$/) {
            ($name, $count) = ($1, $2);
        } else {
            die "Sequence header did not match: $header";
        }
        # these do not wrap!
        chomp $seq;
        #...your loop to read the data goes here
        $add->execute($name, $count, $seq);
    }
}
