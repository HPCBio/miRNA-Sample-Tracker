package miRNA::DB;

use 5.016;
use strict;
use warnings;
#use parent qw(Exporter);
#our @EXPORT_OK=qw();
use List::Compare;
use Scalar::Util qw(blessed);

use DBI;
use Text::CSV_XS;
use Carp;

sub new {
    my ($class, $file) = @_;
    $class = ref $class if ref $class && blessed $class;
    my $self = bless {}, $class;
    $self->_init($file);
}

sub _init {
    my ($self, $dbfile) = @_;
    if (!-e $dbfile) {
        Carp::croak "Database $dbfile doesn't exist";
    }

    my $dbh = DBI->connect("dbi:SQLite:$dbfile","","",{RaiseError =>1})
        or Carp::croak "Couldn't connect to database: " . DBI->errstr;
    $dbh->do("PRAGMA foreign_keys = ON");

    $dbh->do('PRAGMA synchronous = 0');      # Non transaction safe!!!
    $dbh->do('PRAGMA cache_size = 8000000'); # 2 GB dynamic cache increase
                                             # makes index creation faster
    $self->{dbh} = $dbh;
    $self;
}

sub get_sample_names {
    my ($self) = @_;
    $self->{dbh}->selectcol_arrayref(<<SQL);
SELECT name FROM sample
SQL
}

sub get_sample_map {
    my ($self) = @_;
    $self->{dbh}->selectall_hashref(<<SQL, 'sample_id');
SELECT sample_id, name FROM sample
SQL
}

sub get_sample_list {
    my ($self, $list) = @_;
    my @ids;
    {
        local $/ = undef;
        open(my $fh,'<', $list) or die $!;
        @ids = grep { $_ } split(/\n/, <$fh>);
    }
    # run sanity check; every ID must be present and accounted for
    my $sample_names = $self->get_sample_names();
    my $lc = List::Compare->new(\@ids, $sample_names);
    
    if (!$lc->is_LequivalentR()) {
        die "Values don't match sample list in database:input:".
            join(',', $lc->get_unique()).
            "\ndatabase:".join(',', $lc->get_complement());
    }
    
    \@ids;
}

# this will soon be added as a separate table (or tables), but all text for now
sub get_sample_info {
    my ($self, $file) = @_;
    open(my $fh, '<', $file) or die $!;
    
    my $csv = Text::CSV_XS->new({sep_char => "\t"});
    
    my @cols = @{$csv->getline($fh)};
    #my $row = {};
    #$csv->bind_columns (\@{$row}{@cols});
    my %data;
    while (my $d = $csv->getline($fh)) {
        @{$data{$d->[0]}}{@cols[1..$#cols]} = @{$d}[1..$#{$d}];
    }
    \%data, \@cols;
}

#sub get_dbh {
#    my $dbfile = shift;
#    if (!-e $dbfile) {
#        Carp::croak "Database $dbfile doesn't exist";
#    }
#
#    my $dbh = DBI->connect("dbi:SQLite:$dbfile","","",{RaiseError =>1})
#        or Carp::croak "Couldn't connect to database: " . DBI->errstr;
#    $dbh->do("PRAGMA foreign_keys = ON");
#
#    $dbh->do('PRAGMA synchronous = 0');      # Non transaction safe!!!
#    $dbh->do('PRAGMA cache_size = 8000000'); # 2 GB dynamic cache increase
#                                             # makes index creation faster
#    $dbh;
#}

sub get_dbid {
    my ($self, $ref) = @_;
    Carp::croak "No ref given" if !$ref;
    my $id = $self->{dbh}->selectcol_arrayref(<<SQL);
SELECT db_id FROM db_info
WHERE
    db_info.name="$ref"
SQL
    Carp::croak "No ID returned for $ref" if !defined($id) || @$id == 0;
    @$id[0];
}

sub get_db_names {
    my ($self) = @_;
    my $dbnames = $self->{dbh}->selectcol_arrayref(<<SQL);
SELECT name FROM db_info
SQL
    $dbnames;
}

# TODO: this really should go into the database, but it requires a decent amount of de-norm.

sub get_seq_ids {
    my ($self, $reference, $cutoff) = @_;

    # get sequences from specified database based on minimal count
    
    my $dbid = $self->get_dbid($reference);
    
    # TODO: optimize into a column call
    my $mapped_cts = $self->{dbh}->prepare(<<SQL);
SELECT s.srna_id
    FROM srna AS s
    JOIN sample2srna AS sm2s ON (s.srna_id=sm2s.srna_id)
    JOIN srna2db AS s2d ON (s.srna_id=s2d.srna_id)
    WHERE
        s2d.db_id = ?
    AND
        sm2s.count >= ?
    GROUP BY
        s.srna_id
SQL

    my ($id, $count, $seq);
    
    $mapped_cts->bind_columns(\$id);

    $mapped_cts->execute($dbid, $cutoff) or die $mapped_cts->errstr;
    
    my @ids;
    while ($mapped_cts->fetch()) {
        push @ids, $id
    }
    \@ids;
}

sub get_seq_ids_by_total_count {
    my ($self, $reference, $tc) = @_;

    # We ignore $reference here, we only care about total count of samples,
    # which is stored in srna.counts. This is kept for consistency with the
    # above
    my $dbid = $self->get_dbid($reference);
    
    # TODO: optimize into a column call
    my $mapped_cts = $self->{dbh}->prepare(<<SQL);
SELECT s.srna_id
    FROM srna AS s
    JOIN srna2db AS s2d ON (s.srna_id=s2d.srna_id)
    WHERE
        s2d.db_id = ?
    AND
        s.count >= ?
SQL

    my ($id);
    
    $mapped_cts->bind_columns(\$id);

    $mapped_cts->execute($dbid, $tc) or die $mapped_cts->errstr;
    
    my @ids;
    while ($mapped_cts->fetch()) {
        push @ids, $id
    }
    \@ids;
}


# get all counts of a particular sequence across samples
sub get_seq_counts {
    my ($self, $seqid, $reference) = @_;
    my $dbid = $self->get_dbid($reference);
    
    $seqid =~ s/srna//;
    
    # TODO: optimize into a column call
    my $mapped_cts = $self->{dbh}->prepare_cached(<<SQL);
SELECT sm.name, sm2s.count
    FROM sample2srna as sm2s
    JOIN sample as sm ON (sm2s.sample_id=sm.sample_id)
    WHERE
        sm2s.srna_id=?
SQL

    $mapped_cts->execute($seqid) or die $mapped_cts->errstr;
    
    my %cts;
    for my $ct (@{$mapped_cts->fetchall_arrayref()}) {
        $cts{ $ct->[0] } = $ct->[1]
    }
    \%cts;
}

sub get_seq_info {
    my ($self, $seqid, $reference) = @_;
    my $dbid = $self->get_dbid($reference);
    
    $seqid =~ s/srna//;
    
    # TODO: optimize into a column call
    my $mapped_cts = $self->{dbh}->prepare_cached(<<SQL);
SELECT s.srna_id, s.sequence, s2d.hits
    FROM srna AS s
    JOIN srna2db as s2d ON (s.srna_id=s2d.srna_id)
    WHERE
        s.srna_id=?
    AND
        s2d.db_id=?
SQL

    $mapped_cts->execute($seqid, $dbid) or die $mapped_cts->errstr;
    
    my @ids;
    ${$mapped_cts->fetchall_hashref('srna_id')}{$seqid};
}

sub get_db_iterator {
    my ($self, $total_sum, $code) = @_;
    
    #my $dbh = !(ref($db)) ? get_dbh($db) : $db;
    
    $total_sum //= 0;
    
    my $mapped_cts = $self->{dbh}->prepare_cached(<<SQL);
SELECT DISTINCT(s.srna_id), s.sequence, s.count, di.name, s2d.hits
    FROM srna AS s
    LEFT JOIN srna2db as s2d ON (s.srna_id=s2d.srna_id)
    LEFT JOIN db_info as di ON (s2d.db_id=di.db_id)
    WHERE s.count >= ?
SQL

    $mapped_cts->execute($total_sum) or die $mapped_cts->errstr;
    
    my $it = $code ?
        sub { $code->($mapped_cts) } :
        sub { $mapped_cts->fetchrow_hashref; };
    
    return $it;
}

1;

