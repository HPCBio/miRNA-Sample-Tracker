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

    $0 -d <DB_NAME> -c <CUTOFF>

    This script will generate a count table (# times a particular small RNA is
    present in a specific sample).  This can also take a minimal
    cutoff, which represents the minimal count found in any one sample.

USAGE

{

    my $dbfile;
    my $list;
    my $outfile = 'counts.txt';
    my $cutoff = 1;

    GetOptions('db=s'       => \$dbfile,
               'cutoff=i'   => \$cutoff);

    if (!-e $dbfile) {
        die "Database $dbfile doesn't exist";
    }

    my $dbh = DBI->connect("dbi:SQLite:$dbfile","","",{RaiseError =>1})
        or die "Couldn't connect to database: " . DBI->errstr;
    $dbh->do("PRAGMA foreign_keys = ON");

    $dbh->do('PRAGMA synchronous = 0');      # Non transaction safe!!!
    $dbh->do('PRAGMA cache_size = 8000000'); # 2 GB dynamic cache increase
                                             # makes index creation faster

    counts($dbh, $cutoff);

}

sub counts {
    my ($dbh, $cutoff) = @_;

    # need this to get all possible samples (so columns are known in advance)
    
    my $sample_names = get_sample_names($dbh);

    my $get_srnaids = $dbh->prepare(<<SQL);
SELECT srna_id, count FROM srna
WHERE
    count >= ?
SQL

    my $get_cts = $dbh->prepare(<<SQL);
SELECT s.name, s2r.count
    FROM sample2srna s2r
    JOIN sample s ON (s2r.sample_id=s.sample_id)
    WHERE srna_id=?
SQL
    my ($srnaid, $srnaid_ct);

    $get_srnaids->bind_columns(\$srnaid, \$srnaid_ct);
    
    $get_srnaids->execute($cutoff) or die $get_srnaids->errstr;

    my $ct = 0;
    say join("\t", "small_RNA_ID", @{ $sample_names }, 'Sum_Total');
    while ($get_srnaids->fetch()) {
        $get_cts->execute($srnaid) or die $get_cts->errstr;
        my $d = $get_cts->fetchall_hashref('name');
        my @vals = map { exists $d->{$_} ? $d->{$_}->{count} : 0 } @{ $sample_names };
        say join("\t", "srna$srnaid", @vals, (sum @vals));
    }
}

sub get_sample_names {
    my ($dbh) = @_;
    $dbh->selectcol_arrayref(<<SQL);
SELECT name FROM sample
SQL
}
