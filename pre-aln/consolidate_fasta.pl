#!/usr/bin/perl -w
use 5.010;
use strict;
use warnings;
use Time::Piece;
use File::Temp;
use File::Spec;
use Getopt::Long;
use Data::Dumper;

use FindBin;
use lib "$FindBin::Bin/../lib";
use lib $ENV{GIT_REPO}."/hpcbio/smalheiser_scripts/v2/lib";

use miRNA_DB;

my ($db, $new_seq);

GetOptions( 
    'db:s'      => \$db,
    'new:s'     => \$new_seq
);

my $dbh;

if (!$db) {
    die "Need to pass in to -db the old SQLite database"
} else {
    $dbh = miRNA_DB::get_dbh($db);
}

my $it;

if (!$new_seq) {
    die "Need to pass in to -db the old SQLite database"
} else {
    $it = fasta_iterator($new_seq);
}


my $highest = miRNA_DB::last_seqid($dbh);

# this tries to find a sane start point for any new sequence IDs, basically just
# incrementing the first digit and buffering with the appropriate # of 0's
(my $next_start = $highest) =~ s/^(\d)(\d+)/ ($1 + 1).('0' x length($2)) /e;

my $found = my $uniq = 0;

#say "$highest -> $next_start";
#
#exit;

(my $newly_collapsed = $new_seq) =~ s/\.\w+?$/.consolidated.fasta/;

open(my $new_fh, '>', $newly_collapsed) or die $!;

my $ct = 0;
while (my $fasta = $it->()) {
    last if !$fasta->{seq};
    my $id = miRNA_DB::get_seqid($dbh, $fasta->{seq});
    if ($id) {
        $found++;
    } else {
        $uniq++;
        $id = $next_start++;
    }
    say $new_fh ">srna$id ".$fasta->{count}."\n".$fasta->{seq};
    #last if $ct++ == 100000;
}

say "Found: $found\nUnique:$uniq";

sub fasta_iterator {
    my $file = shift;
    open(my $fh, '<', $file) or die $!;
    return sub {
        local $/ = "\n>";
        while (my $chunk = <$fh>) {
            $chunk =~ tr/>//d;
            my ($header, $seq) = split("\n", $chunk, 2);
            $seq =~ s/\s+//g;
            my ($orig_id, $count);
            if ($header =~ /^(\S+)\s(\d+)$/) {
                ($orig_id, $count) = ($1, $2);
            } else {
                die $header;
            }
            return {
                'orig_id'   => $orig_id,
                'count'     => $count,
                'seq'       => $seq
            };
        }
    }
}