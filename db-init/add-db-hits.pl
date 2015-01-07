#!/usr/bin/env perl
use 5.010;
use strict;
use warnings;
use autodie;

use lib $ENV{GIT_REPO}."/hpcbio/smalheiser_scripts/v2/lib";
use FindBin;
use lib "$FindBin::Bin/../lib";

use File::Spec;
use File::Basename;
use Getopt::Long;
use DBI;

my $USAGE = <<USAGE;

    $0 -d <DB_NAME> -r <REFERENCE_NAME> -b <BAM_FILE>

    This script will read a BAM file and add simple information to the SQLite3
    database regarding the name of the reference database and # of hits a
    specific small RNA read had against the reference database. No other
    information is stored in the database at this time, as any location-specific
    information can be directly accessed via the BAM file.

USAGE

{

    my $bam;
    my $dbfile;
    my $ref;

    GetOptions('db=s'       => \$dbfile,
               'bam=s'      => \$bam,
               'ref=s'      => \$ref);

    die $USAGE unless (defined $ref) && (defined $bam && -e $bam);

    if (!-e $dbfile) {
        die "Database $dbfile doesn't exist";
    }

    my $dbh = DBI->connect("dbi:SQLite:$dbfile","","",{RaiseError =>1})
        or die "Couldn't connect to database: " . DBI->errstr;
    $dbh->do("PRAGMA foreign_keys = ON");

    $dbh->do('PRAGMA synchronous = 0');      # Non transaction safe!!!
    $dbh->do('PRAGMA cache_size = 8000000'); # 2 GB dynamic cache increase
                                             # makes index creation faster

    say STDERR "Starting BAM parse";
    $dbh->do("BEGIN");
    load_hits($dbh, $bam, $ref);
    $dbh->do("COMMIT");
    say STDERR "Finished BAM parse";

    # TODO: add a db-specific index for performance?

    #$dbh->do("BEGIN");
    #$dbh->do ("CREATE INDEX srna2db_db_idx on srna2db( db_id )");
    #$dbh->do("COMMIT");

}

sub load_hits {
    my ($dbh, $bam, $ref) = @_;

    my $db_id = add_sample($dbh, $ref);

    my $add_hits = $dbh->prepare(<<SQL);
INSERT INTO srna2db (srna_id, db_id, hits)
    VALUES(?, ?, ?)
SQL

    open(my $BAM_STREAM, "samtools view -F 4 $bam | ");

    my %seen;
    my $ct = 0;
    say STDERR "Loading $bam into hits table";
    while ( <$BAM_STREAM> ) {
        my ($read, @rest) = split("\t", $_, 12);
        next if $seen{$read}++;
        my $nh = 1;
        if ($rest[-1] =~ /NH:i:(\d+)/ ) {
            $nh = $1;
        }
        my $id;
        if ($read =~ /srna(\d+)/) {
            $id = $1;
        } else {
            die "Problem with line $.: $_";
        }
        $add_hits->execute($id, $db_id, $nh) or die $add_hits->errstr;
        $ct++;
        print "Loaded $ct records\n" if ! $ct % 100000;
    }
}

sub add_sample {
    my ($dbh, $nm) = @_;
    say STDERR "Loading $nm into database table";
    my $add_sample = $dbh->prepare(<<SQL);
INSERT OR IGNORE INTO db_info ( name )
    VALUES(?)
SQL
    $add_sample->execute($nm);
    my $id = $dbh->last_insert_id('','','','');
    die "No ID returned, maybe $nm is already present in the database?" unless $id;

    $id;
}
