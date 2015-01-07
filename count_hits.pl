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
use List::Util qw(sum);

my $USAGE = <<USAGE;

    $0 -d <DB_NAME> -l <ID LIST>

    This script will generate a database hit table (# times a particular small
    RNA has hits in a particular database).

TODO: Take an optional ID list (file with a list of small RNA IDs, separated by newlines).

USAGE

{

    my $dbfile;
    my $list;
    my $outfile = 'counts.txt';

    GetOptions('db=s'       => \$dbfile,
               'list=s'     => \$list);

    if (!-e $dbfile) {
        die "Database $dbfile doesn't exist";
    }

    my $dbh = DBI->connect("dbi:SQLite:$dbfile","","",{RaiseError =>1})
        or die "Couldn't connect to database: " . DBI->errstr;
    $dbh->do("PRAGMA foreign_keys = ON");

    $dbh->do('PRAGMA synchronous = 0');      # Non transaction safe!!!
    $dbh->do('PRAGMA cache_size = 8000000'); # 2 GB dynamic cache increase
                                             # makes index creation faster

    counts($dbh);

}

sub counts {
    my ($dbh) = @_;

    # need this to get all possible samples (so columns are known in advance)

    my $db_names = get_sample_names($dbh);

    my $get_srnaids = $dbh->prepare(<<SQL);
SELECT srna_id, count FROM srna
SQL

    my $get_cts = $dbh->prepare(<<SQL);
SELECT s2d.srna_id, d.name, s2d.hits
    FROM srna2db s2d
    JOIN db_info d ON (s2d.db_id=d.db_id)
    WHERE srna_id=?
SQL
    my ($srnaid, $srnaid_ct);

    $get_srnaids->bind_columns(\$srnaid, \$srnaid_ct);

    $get_srnaids->execute() or die $get_srnaids->errstr;

    my $ct = 0;
    say join("\t", "small_RNA_ID", 'Count', @{ $db_names });
    while ($get_srnaids->fetch()) {
        $get_cts->execute($srnaid) or die $get_cts->errstr;
        my $d = $get_cts->fetchall_hashref('name');
        my @vals = map { exists $d->{$_} ? $d->{$_}->{hits} : 0 } @{ $db_names };
        say join("\t", "srna$srnaid", $srnaid_ct, @vals);
    }
}

sub get_sample_names {
    my ($dbh) = @_;
    $dbh->selectcol_arrayref(<<SQL);
SELECT name FROM db_info
SQL
}
