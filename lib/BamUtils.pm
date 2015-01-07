package BamUtils;

use 5.016;
use strict;
use warnings;
no if $] >= 5.017011, warnings => 'experimental::smartmatch';
use feature "switch";
use Data::Dumper;
#use Bio::DB::Sam;
use parent qw(Exporter);

our @EXPORT_OK=qw();

{
    
    my %cache;
sub flag_stats {
    my ($flag) = @_;
    my %stats;
    
    if (exists $cache{$flag}) {
        return $cache{$flag};
    }
    
    for ($flag) {
        if ($_ & 0x1) {
            $stats{multiple_segs}++;
        }
        if ($_ & 0x2) {
            $stats{properly_aligned}++;
        }
        if ($_ & 0x4) {
            $stats{unmapped}++;
        }
        if ($_ & 0x8) {
            $stats{next_seg_unmapped}++;
        }
        if ($_ & 0x10) {
            $stats{revcomp}++;
        }
        if ($_ & 0x20) {
            $stats{next_seg_revcomp}++;
        }
        if ($_ & 0x40) {
            $stats{first_seg}++;
        }
        if ($_ & 0x80) {
            $stats{last_seg}++;
        }
        if ($_ & 0x100) {
            $stats{secondary_aln}++;
        }
        if ($_ & 0x200) {
            $stats{no_pass_QC}++;            
        }
        if ($_ & 0x400) {
            $stats{'PCR/optical_dup'}++;
        }
        if ($_ & 0x800) {
            $stats{'supl_aln'}++;
        }
    }
    
    $cache{$flag} = \%stats;
    
    return \%stats;
}
}

{
    
    my %str_cache;
sub flag_strand {
    my ($flag) = @_;
    my $strand = '+';
    
    if (exists $str_cache{$flag}) {
        return $str_cache{$flag};
    }
    
    for ($flag) {
        if ($_ & 0x10) {
            $strand = '-';
        }
    }
    
    $str_cache{$flag} = $strand;
    
    return $strand;
}
}

sub cigar_stats {
    my ($cigar) = @_;
    my $type;
    if ($cigar =~ /^\d+M$/) {
        $type = 'EXACT';
    } elsif ($cigar =~ /^(\d+[HS])?\d+M(\d+[HS])?$/) {
        $type = 'CLIPPED';
    } else {
        $type = 'MISMATCH/INDEL';
    }
    $type;
}

sub cigar_len {
    my ($cigar) = @_;
    my $len = 0;
    while ($cigar =~ /(\d+)[MDN=XP]/g) {
        $len += $1;
    }
    $len;
}

sub tag_info {
    my ($tags) = @_;
    my %ti; 
    while ($tags =~ /((?:Z|N)[NHL]):\w:(\S+)/g) {
        for ($1) {
            when ('ZH') {
                $ti{hpscore} = $2;
            }
            when ('ZN') {
                $ti{alns} = $2;
            }
            when ('NH') {
                $ti{alns} = $2;
            }
            when ('ZL') {
                $ti{hploc} = $2;
            }
        }
    }
    $ti{alns} //= 1;
    \%ti;
}

# using a samtools pipe, get a mapping of ref seqs and locations to the query
# name from a BAM file and a list of wanted IDs, returns a mapping to this for
# downstream use. Works best with a query-name sorted BAM and a sum cutoff (I'm
# sure this would tank with no cutoff)

sub get_ref_map {
    my ($bam, $ids) = @_;  # BAM file and ID map
    
    die "ID map must be an array ref" if $ids && ref $ids ne 'ARRAY';
    my %in_ids = map { "srna".$_ => $_ } @$ids;
    
    my %map;
    open(my $samfh, "samtools view -F 4 $bam | ") or die $!;
    
    my $ct = 0;
    while (<$samfh>) {
        chomp;
        my $d = [split("\t", $_, 12)];
        
        next unless exists $in_ids{$d->[0]};
        
        # simple hash of the IDs vs reference hits (note there is a full
        # location option below commented out)
        $map{ $in_ids{ $d->[0] } }{ $d->[2] }++;
        #push @{ $map{ $in_ids{$d->[0]} } }, sprintf("%s:%d-%d(%d)", $d->[2], $d->[3], $d->[3] + cigar_len($d->[5]) - 1, ($d->[1] & 0x10) ? -1 : 1 );
    }
    \%map;
}

# these need to be replaced with Bio-Samtools variants, but that doesnt work
# with latest version yet. We punt for now.
sub get_bam_iterator {
    my ($bam, $code) = @_;

    check_query_sorted($bam);
    
    open(my $samfh, "samtools view -F 4 $bam | ") or die $!;
    
    my $cb = sub {
        my $l = <$samfh>;
        return unless $l;
        process_bam_line($l);
    };
    
    # allowing a piece of code to process the iterator
    # enables chunking code together
    my $it = $code ?
        sub { $code->($cb) } :
        sub { $cb->(); };
    
    return $it;
}

sub check_query_sorted {
    my ($bam) = @_;
    
    open(my $samfh, "samtools view -H $bam | ") or die $!;
    my $so;
    while (<$samfh>) {
        if (/\s+SO:\s*(\S+)/) {
            $so = $1;
            last;
        }
    }
    
    if ( !$so ) {
        die "BAM file does not have sort order"
    }
    
    if ( $so !~ /query/i) {
        die "BAM file is not query name-sorted, SO: $so"
    }
    
    1;
}

sub process_bam_line {
    my $l = shift;
    my @data = split("\t", $l, 12);
    my %tags = map {
        my ($nm, $type, $data) = split(":", $_, 3);
        (
            $nm => { type  => $type,
                     data  => $data }
        )
        } split("\t", $data[11]);
    $data[11] = \%tags;
    \@data;
}

1;