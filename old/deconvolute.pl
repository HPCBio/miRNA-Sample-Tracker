#!/usr/bin/env perl
use 5.010;
use strict;
use warnings;
use autodie;
use File::Spec;
use Getopt::Long;
use DBI;

my $USAGE = <<USAGE;

    $0 -f <sample.compressed.fasta> -db <DB_NAME> --dir <DIR> -c

This reads a FASTA from a collapsed sample sequence file and a database
(created using fastq2tbl.pl) and prints a simple tab-delimited list of read ID
and counts for that sample. Match is based upon the sequence read.

If given a directory this will simply iterate through all the files and
(sequentially) determine counts for each.  It ridiculously NON-parallel at this
time.

TODO: optimize this into a table for the database?

USAGE


{

    my ($fasta, $dir);
    my $dbfile = 'srna.sqlite3';
    my $create = 0;

    GetOptions('fasta=s'    => \$fasta,
               'db=s'       => \$dbfile,
               'dir=s'      => \$dir,
               'create'     => \$create);

    die unless (defined $fasta && -e $fasta) || (defined $dir && -e $dir);

    # TODO: allow other extensions?
    my @fastas = $dir ? glob(File::Spec->catfile($dir, '*.derep.fasta')) : $fasta;

    if (!-e $dbfile) {
        die "Database $dbfile doesn't exist";
    }

    my $dbh = DBI->connect("dbi:SQLite:$dbfile","","",{RaiseError =>1})
        or die "Couldn't connect to database: " . DBI->errstr;

    for my $f (@fastas) {
        get_data($dbh, $f);
    }

}

sub get_data {
    my ($dbh, $fasta_file) = @_;
    my $sth = $dbh->prepare(<<SQL);
SELECT name
    FROM srna
    WHERE sequence = ?;
SQL
    open(my $FASTA, '<', $fasta_file);
    open(my $COUNTS, '>', "$fasta_file.counts");

    my ($nm);

    $sth->bind_columns(\$nm);
    while (!eof($FASTA)) {
        my $header = <$FASTA>;
        my $seq = <$FASTA>;
        my ($name, $count);
        if ($header =~ /^>(\d+)-(\d+)$/) {
            ($name, $count) = ($1, $2);
        } else {
            die "Sequence header did not match: $header";
        }
        # these do not wrap!
        chomp $seq;
        $sth->execute($seq) or die $sth->errstr;
        $sth->fetch() or die $sth->errstr;
        say $COUNTS join("\t", $nm, $count);
    }
    close $FASTA;
    close $COUNTS;
    1;
}
