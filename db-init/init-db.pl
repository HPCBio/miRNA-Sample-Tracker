#!/usr/bin/env perl
use 5.010;
use strict;
use warnings;
use autodie;

use FindBin;
use lib "$FindBin::Bin/../lib";

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
my $dbfile = 'test.sqlite3';
my $create = 1;

GetOptions('fasta=s'    => \$fasta,
           'db=s'       => \$dbfile,
           'create'     => \$create);

if (-e $dbfile) {
    die "Database $dbfile exists, set --create to overwrite" if !$create;
    unlink $dbfile;
}

my @ORDER = qw(srna sample sample2srna db_info srna2db);

my %SCHEMA = (
    srna    => <<SQL,
    (
    srna_id   INTEGER(10) PRIMARY KEY,
    count     INTEGER(10),     -- counts for this sequence in all samples, may be null
    sequence  VARCHAR(50) NOT NULL
    );
SQL

    sample => <<SQL,
    (
    sample_id  INTEGER PRIMARY KEY AUTOINCREMENT,
    name       VARCHAR(50),
    raw_reads  INTEGER(10)
    );
SQL

    sample2srna => <<SQL,
    (
    srna_id      INTEGER(10) REFERENCES srna,
    sample_id    INTEGER(10) REFERENCES sample,
    count        INTEGER(10),
    PRIMARY KEY (srna_id, sample_id)
    );
SQL

    db_info => <<SQL,
    (
    db_id     INTEGER PRIMARY KEY AUTOINCREMENT,
    name      VARCHAR(50) UNIQUE
    );
SQL

    srna2db => <<SQL,
    (
    srna_id      INTEGER(10) REFERENCES srna,
    db_id        INTEGER(10) REFERENCES db_info,
    hits         INTEGER(10),
    PRIMARY KEY (srna_id, db_id)
    );
SQL
);


my $dbh = DBI->connect("dbi:SQLite:$dbfile","","",{RaiseError =>1})
    or die "Couldn't connect to database: " . DBI->errstr;

$dbh->do("PRAGMA foreign_keys = ON");

if ( $create ) {
    print "Creating new DNA Database\n";
    for my $table (@ORDER) {
        die unless exists $SCHEMA{$table};
        #say "CREATE TABLE $table $SCHEMA{$table}";
        $dbh->do("CREATE TABLE $table $SCHEMA{$table}");
    }
}

if ($fasta && -e $fasta) {

    $dbh->do('PRAGMA synchronous = 0');      # Non transaction safe!!!
    $dbh->do('PRAGMA cache_size = 4000000'); # 2 GB dynamic cache increase
                                             # makes index creation faster

    say STDERR "Loading small RNA data";
    $dbh->do("BEGIN");
    import_srna($dbh, $fasta);
    $dbh->do("COMMIT");

    say STDERR "Indexing small RNA data";
    $dbh->do("BEGIN");
    $dbh->do ("CREATE UNIQUE INDEX srna_seq_idx ON srna (sequence)");
    $dbh->do("COMMIT");

    say STDERR "Finished small RNA data";

}

exit;

}

sub import_srna {
    my ($dbh, $fasta_file) = @_;

    # TODO: pass this statement into the method instead of the dbh
    my $add = $dbh->prepare("INSERT INTO srna ( srna_id, count, sequence ) VALUES(?,?,?)");

    open(my $FASTA, '<', $fasta_file);

    while (!eof($FASTA)) {
        my $header = <$FASTA>;
        my $seq = <$FASTA>;
        my ($name, $count);
        if ($header =~ /^>srna(\d+)\s(\d+)$/) {
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
