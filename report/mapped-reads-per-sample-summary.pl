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

my $USAGE = <<USAGE;

    $0 -d <DB_NAME> -c <CUTOFF>

    This script will generate a summary table of the total number of mapped hits
    per sample, taking into account the sample sequence redundancy. User must
    provide the SQLite3 database. An optional cutoff level can be added
    (default=1, or keep everything)

USAGE

{

    my $dbfile;
    #my $cutoff = 1;
    #my $sum = 0;

    GetOptions('db=s'       => \$dbfile,
               #'cutoff=i'   => \$cutoff,
               #'sum=i'      => \$sum
               );

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
    #my $db_id = get_db_names($dbh, $ref);

    # reads mapping to the database
    my $mapped_cts = $dbh->prepare(<<SQL);
SELECT sm.name, sm.sample_id, SUM(sm2s.count) AS Total
    FROM sample AS sm
    JOIN sample2srna AS sm2s ON (sm.sample_id=sm2s.sample_id)
    WHERE
        sm2s.srna_id IN (
            SELECT s2d.srna_id
                FROM srna2db AS s2d
                WHERE s2d.db_id=?
        )
    GROUP BY sm.sample_id
SQL

    my $mapped_any_cts = $dbh->prepare(<<SQL);
SELECT sm.name, sm.sample_id, SUM(sm2s.count) AS Total
    FROM sample AS sm
    JOIN sample2srna AS sm2s ON (sm.sample_id=sm2s.sample_id)
        sm2s.srna_id IN (
            SELECT s2d.srna_id
                FROM srna2db AS s2d
        )
    GROUP BY sm.sample_id
SQL

    my $all_cts = $dbh->prepare(<<SQL);
SELECT sm.name, sm.sample_id, SUM(sm2s.count) AS Total
    FROM sample AS sm
    JOIN sample2srna AS sm2s ON (sm.sample_id=sm2s.sample_id)
    GROUP BY sm.sample_id
SQL

    open(my $outfh, '>', "mapping.txt") or die $!;

    my %summary;
    
    my @refs = get_db_names($dbh);
    
    my ($nm, $sid, $total);
    for my $r (@refs) {
        my ($db_id, $refname) = @$r;
        $mapped_cts->bind_columns(\$nm, \$sid, \$total);
        
        $mapped_cts->execute($db_id) or die $mapped_cts->errstr;
        while ($mapped_cts->fetch()) {
            @{ $summary{$sid}{$refname} }{ qw(Mapped) } = ($total);
        }
        say STDERR "Finished $refname";
    }
    
    my @refnames = map { $_->[1] } @refs;
    
    $mapped_any_cts->bind_columns(\$nm, \$sid, \$total);
        
    $mapped_any_cts->execute() or die $mapped_any_cts->errstr;
    
    while ($mapped_any_cts->fetch()) {
        @{ $summary{$sid}{any_mapped} }{ qw(Mapped) } = ($total);
    }
    
    say STDERR "Finished any_mapped";

    say $outfh join("\t", "Name", "Raw_Cutoff", map {
        my $id  = $_;
        join("\t", map {"$id $_"} ('Mapped', 'Percent Mapped'));
        } (@refnames, 'any_mapped'));
    
    $all_cts->bind_columns(\$nm, \$sid, \$total);
    $all_cts->execute() or die $all_cts->errstr;
    
    while ($all_cts->fetch()) {
        #say STDERR join("\t", $nm, $sid, $total);
        my @vals = ($nm, $total);
        
        say $outfh join("\t",
                        $nm,
                        $total,
                        map { my $id = $_;
                             join("\t",
                                  $summary{$sid}{$id}{Mapped},
                                  sprintf("%.2f", 100 * ($summary{$sid}{$id}{Mapped}/$total))
                                  )
                             } (@refnames, 'any_mapped'));
    }
}

sub get_db_names {
    my ($dbh, $ref) = @_;
    my $id = $dbh->selectall_arrayref(<<SQL);
SELECT db_id, name FROM db_info
SQL
    @$id;
}
