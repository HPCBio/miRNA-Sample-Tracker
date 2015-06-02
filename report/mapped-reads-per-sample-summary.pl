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

    $0 -d <DB_NAME> -format <tab/textile/md>

    This script will generate a summary table of the total number of mapped hits
    per sample, taking into account the sample sequence redundancy. User must
    provide the SQLite3 database.

USAGE

{

    my %VALID = (
        tab         => {
                        'heading'   => \&to_tab,
                        'header'    => \&to_tab,
                        'row'       => \&to_tab,
                       },
        md          => {
                        'heading'   => \&md_heading,
                        'header'    => \&md_tbl_header,
                        'row'       => \&generic_row,
                       },
        textile     => {
                        'heading'   => \&textile_heading,
                        'header'    => \&textile_tbl_header,
                        'row'       => \&generic_row,
                       }
    );

    my $dbfile;
    my $format = 'tab';
    my $prefix;

    GetOptions('db=s'       => \$dbfile,
               'format:s'   => \$format,
               'prefix:s'   => \$prefix
               #'cutoff=i'   => \$cutoff,
               #'sum=i'      => \$sum
               );

    if (!-e $dbfile) {
        die "Database $dbfile doesn't exist";
    }

    if (! exists $VALID{ $format }) {
        die "Format $format not supported\n\n$USAGE";
    }
    
    my $callbacks = $VALID{ $format };
    
    my $dbh = DBI->connect("dbi:SQLite:$dbfile","","",{RaiseError =>1})
        or die "Couldn't connect to database: " . DBI->errstr;
    $dbh->do("PRAGMA foreign_keys = ON");

    $dbh->do('PRAGMA synchronous = 0');      # Non transaction safe!!!
    $dbh->do('PRAGMA cache_size = 8000000'); # 2 GB dynamic cache increase
                                             # makes index creation faster
    
    counts($dbh, $callbacks);
}

sub counts {
    my ($dbh, $fmt, $outfh) = @_;

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
    WHERE
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

    $fmt->{heading}('Mapping summary');
    
    $fmt->{header}(
        ["Name",
         "Total Reads",
         map {
        my $id  = $_;
        map {"$id $_"} ('Mapped', 'Percent Mapped');
        } (@refnames, 'any_mapped')]);
    
    $all_cts->bind_columns(\$nm, \$sid, \$total);
    $all_cts->execute() or die $all_cts->errstr;
    
    while ($all_cts->fetch()) {
        #say STDERR join("\t", $nm, $sid, $total);
        my @vals = ($nm, $total);
        
        $fmt->{row}([
                        $nm,
                        $total,
                        map { my $id = $_;
                 $summary{$sid}{$id}{Mapped},
                 sprintf("%.2f", 100 * ($summary{$sid}{$id}{Mapped}/$total)
                 )
            } (@refnames, 'any_mapped')]);
    }
}

sub get_db_names {
    my ($dbh, $ref) = @_;
    my $id = $dbh->selectall_arrayref(<<SQL);
SELECT db_id, name FROM db_info
SQL
    @$id;
}

sub to_tab {
    say join("\t", ref $_[0] eq 'ARRAY' ? @{$_[0]} : $_[0])
}

sub md_heading {
    my ($str) = @_;
    say "## $str";
}

sub md_tbl_header {
    my ($colnames) = @_;
    say '|'.join("|", @$colnames).'|';
    say '|'.join("|",
        map {
            ':'.
            ('-' x (length($_) - 1))
            }
        @$colnames ).'|'
}

sub textile_heading {
    my ($str) = @_;
    say "h2. $str\n";
}

sub textile_tbl_header {
    my ($colnames) = @_;
    #say '|^.';
    say '|_. '.join("|_. ", @$colnames).'|';
    #say '|-.';
}

sub generic_row {
    my ($data) = @_;
    say '|'.join("|", @$data).'|';
}