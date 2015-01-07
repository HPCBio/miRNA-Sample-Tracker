#!/usr/bin/env perl
use 5.014;
use strict;
use warnings;
use autodie;
use Data::Dumper;
use Getopt::Long;
use Text::CSV_XS;
use lib $ENV{GIT_REPO}."/hpcbio/smalheiser_scripts/v2/lib";
use miRNA_DB;
use FeatureUtils;
use SampleUtils;
use File::Basename;
use Carp;

my %VALID_FORMATS = map {$_ => 1} qw(gff bed);

my @VALID_SOURCES = qw(wgRNA tRNAs lincRNAsTranscripts refseq_canon);

{
    
my $bam;
my $features;
my $format = 'bed';
my $db;
my $cutoff = 0;
my $ref = 'hg19_spiked';
my $kgXref = '/home/mirrors/igenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/kgXref.txt';
my $max = 0;
my $min = 1;
my $sum = 0;
my $addcts = 0;
my $meta_file;
my $sample_order_file;

my $USAGE = <<USAGE;

    perl $0 -b <BAM> -f <ANNOTATION.bed/gff> -t <bed/gff> -o my_table.txt -d <DATABASE> -c <CUTOFF> -r <REFERENCE>

Requires -b, -f, -d, -k

Also requires modules 'perl/5.16.1_hpcbio', 'samtools', and 'bedtools2'

Works best with a queryname-sorted BAM file, which will sort the srna IDs in order

Ex: perl smalheiser_scripts/v2/smal/intersect-to-tables.pl \
        -b hg19_spike_ALL/all.name.bam \
        -f ../../annotation/all.gtf \
        -d test.sqlite3 \
        -r hg19_spiked \
        -max 10 \
        -kg ../../annotation/kgXref.txt \
        -sum 500 \
        -a \
        -l sample_list.txt \
        -meta sample_info2.txt

USAGE

GetOptions(
    'bam=s'         => \$bam,
    'features=s'    => \$features,
    'type=s'        => \$format,
    'sum_co=s'      => \$sum,      # total count filter (collides with cutoff)
    'database=s'    => \$db,
    'cutoff=i'      => \$cutoff,  # min count for srna in any one sample
    'reference=s'   => \$ref,
    'kgXref=s'      => \$kgXref,
    'min=i'         => \$min,
    'max=i'         => \$max,     # max # hits in reference; may add max # alns reported as well
    'addcts'        => \$addcts,
    'meta=s'        => \$meta_file,
    'list=s'        => \$sample_order_file
);

#my $source = File::Basename::basename($features);

die $USAGE if !defined($bam) || !defined($features) || !defined($db) || !defined($kgXref) ;

die $USAGE if $sum && $cutoff; # can't define both

my $dbxrefs = FeatureUtils::kgXref($kgXref);

my $csv = Text::CSV_XS->new({
                             sep_char  => "\t",
                             allow_loose_quotes => 1
                             });

my $dbh = miRNA_DB::get_dbh($db);

my $sl = [];

my ($meta, $metacols);

if ($addcts) {
    if ($sample_order_file) {
        $sl = miRNA_DB::get_sample_list($dbh, $sample_order_file);
    } else {
        $sl = miRNA_DB::get_sample_names($dbh);
    }
    
    if ($meta_file) {
        ($meta, $metacols) = SampleUtils::get_sample_info($meta_file);
    }
}

my %cutoff_ids;

if ($sum) {
    %cutoff_ids = map {$_ => 1} @{miRNA_DB::get_seq_ids_by_total_count($dbh, $ref, $sum)};
} else { 
    %cutoff_ids = map {$_ => 1} @{miRNA_DB::get_seq_ids($dbh, $ref, $cutoff)};
}

# samtools view -bh final_filtered.bam | intersectBed -abam - -b bed/trna.bed -wa -wb -bed | sortBed | cut -f1-4,6,14-16,18 > final_filtered.counts.trna.txt
#open(my $in, "samtools view -F 4 -bh $bam | intersectBed -abam - -b $features -loj -wb -bed |");

open(my $in, "samtools view -F 4 -bh $bam | bamToBed -cigar -i stdin | intersectBed -a stdin -b $features -loj -wb |");

## Input:
##                            0      1          2      3    4      5       6 
## a -> BED-ish: col 0-11  chrom chromStart chromEnd qname score strand  cigar
#
##                           7     8      9       10        11      12    13     14     15
## b -> GTF: col 12-20     chrom source stype chromStart chromEnd score strand frame attributes
#
## Output:
## READ SEQ LENGTH NUM_ALNS REF REF_START REF_END REF_STRAND MATCH_TYPE
## SOURCE_FEATURES FEATURE_NAME GENE_SYMBOL FEATURE_START FEATURE_END
## FEATURE_STRAND DESCRIPTION
#
## Requires:
## BAM info (READ SEQ LENGTH NUM_ALNS REF REF_START REF_END REF_STRAND MATCH_TYPE),
## database (NUM_ALNS),
## *.gtf (SOURCE_FEATURES FEATURE_NAME GENE_SYMBOL FEATURE_START FEATURE_END FEATURE_STRAND),
## kgXref (DESCRIPTION and other bits)

my @filtered_cols = qw(READ SEQ LENGTH NUM_ALNS);
my @nongenic_cols = (@filtered_cols, 'HIT_LOCATIONS' );
my @genic_cols = (@nongenic_cols, @VALID_SOURCES, 'RefSeq DESCRIPTION');

my $filter = $sum ? "$sum-totalcount" : "$cutoff-cutoff";

open(my $genic_fh, '>', "genic.$filter.$min-$max-alns.txt" ) or die $!;
open(my $nongenic_fh, '>', "nongenic.$filter.$min-$max-alns.txt" ) or die $!;
open(my $filtered_fh, '>', "filtered.$filter.$min-$max-alns.txt" ) or die $!;

if ($meta) {
    say $genic_fh SampleUtils::meta_data_str($meta, $sl, $#genic_cols);
    say $nongenic_fh SampleUtils::meta_data_str($meta, $sl, $#nongenic_cols);
    say $filtered_fh SampleUtils::meta_data_str($meta, $sl, $#filtered_cols);
}

#say join("\t", @filtered_cols, @$sl);
say $genic_fh join("\t", @genic_cols, @$sl);
say $nongenic_fh join("\t", @nongenic_cols, @$sl);
say $filtered_fh join("\t", @filtered_cols, @$sl);

my @cols = qw(bam.chrom bam.chromStart bam.chromEnd bam.qname bam.score bam.strand bam.cigar
gtf.chrom gtf.source gtf.stype gtf.chromStart gtf.chromEnd gtf.score gtf.strand gtf.frame gtf.attributes);

$csv->column_names(@cols);

my $ct = 0;

my %read_buffer;

my $last_details;

# TODO: invert logic, move into an iterator to both open the input stream and
# chunk data via a callback; this can wrap the Text::CSV stuff 
# e.g. 
# my $it = generate_processor($input_file, \&callback, \%args);
# while ( my $d = $it->() ) {  # output... }

say "Total IDs: ".(scalar keys %cutoff_ids);

# TODO: ugly f*$%ing hack, need to fix into an iterator (see above)
FILTER:
while (my $d = $csv->getline_hr($in)) {
    # filter based on cutoff
    
    ($d->{sid} = $d->{'bam.qname'}) =~ s/srna//;
    next FILTER unless exists $cutoff_ids{ $d->{sid} };
    
    # cludge to clean up source name a little more
    $d->{'gtf.source'} =~ s/hg19_//;

    my ($processed, $seq_data, $status);
    
    if ($last_details && $d->{'bam.qname'} ne $last_details->{'bam.qname'}) {
    
        say Dumper $last_details, $d if !keys %read_buffer;
        # TODO: Could (with a little refactoring) allow limiting the # alns to a random pick: from Randal Scwartz:
        # my @result;
        # push @result, splice @list, rand @list, 1 while @result < $n;
        my ($str, $status) = process_data(\%read_buffer, $dbxrefs, $dbh, $ref, $max, $min) if $last_details->{'bam.qname'};
        
        #say $str;
        
        if ($addcts) {
            my $counts = miRNA_DB::get_seq_counts($dbh, $read_buffer{'bam.qname'}, $ref);
            $str .= "\t".join("\t", map { exists $counts->{$_} ? $counts->{$_} : 0 } @$sl);
        }
        
        # TODO: did I mention this is a convoluted mess?
        if ($status == -1) {
            say $filtered_fh $str;
        } else {
            if ( $status == 0 ) { # no features
                say $nongenic_fh $str;
            } else {
                say $genic_fh $str;
            }
        }
        %read_buffer = ();
        delete $cutoff_ids{ $last_details->{'sid'} };
        
        print "IDs left:".scalar(keys %cutoff_ids);
        if (scalar(keys %cutoff_ids) <= 10) {
            say "\t", join(",", keys %cutoff_ids);
        } else {
            print "\n";
        }
        
    }
    last FILTER if keys %cutoff_ids == 1;
    
    $read_buffer{'bam.qname'} = $d->{'bam.qname'};
    push @{ $read_buffer{data}{$d->{'gtf.source'} } }, $d;
    
    $last_details = $d;
}

# fencepost
if (scalar keys %read_buffer) {
    my ($str, $status) = process_data(\%read_buffer, $dbxrefs, $dbh, $ref, $max, $min);
    if (ref $str) {
            say $filtered_fh join("\t",
                                  $last_details->{'bam.qname'},
                                  $str->{sequence},
                                  length($str->{sequence}),
                                  $str->{hits});
    } else {
        if ($status) {
            say $genic_fh $str;
        } else {
            say $nongenic_fh $str;
        }
    }
}

exit;

}

# TODO: ugly f*$%ing hack, need to fix into an iterator
sub process_data {
    my ($buffer, $kgxref, $dbh, $ref, $max, $min) = @_;
    
    Carp::croak "No SeqID ".Dumper($buffer) unless $buffer->{'bam.qname'};
    
    my $seq_data = miRNA_DB::get_seq_info($dbh, $buffer->{'bam.qname'}, $ref);
    
    my $str;
    
    # TODO: optimize these two checks down
    # TODO: note this should be taken care of at the init. db query level (all
    # relevant info is there) filter if NUM_ALNS > $max
    if ($max && ($seq_data->{hits} > $max)) {
        $str = join("\t",
            $buffer->{'bam.qname'},
            $seq_data->{sequence},
            length($seq_data->{sequence}),
            $seq_data->{hits});

        return ($str, -1);
    }
    
    # filter if NUM_ALNS < $min
    if ($min && ($seq_data->{hits} < $min)) {
        $str = join("\t",
            $buffer->{'bam.qname'},
            $seq_data->{sequence},
            length($seq_data->{sequence}),
            $seq_data->{hits});
        
        return ($str, -1);
    }

    # rearrange to group everything; we don't care about context much here, just uniqueness
    my %features;
    for my $src (@VALID_SOURCES, '.') {
        # we treat refseq_canon differently
        if (exists $buffer->{data}{$src} ) {
            for my $hit (@{ $buffer->{data}{$src} } ) {
                my $loc = sprintf("%s:%d-%d(%s)", $hit->{'bam.chrom'}, $hit->{'bam.chromStart'}, $hit->{'bam.chromEnd'}, $hit->{'bam.strand'});
                $features{'location'}{$loc}++;
                if ($src ne '.') {
                    my $atts = FeatureUtils::gtf_attributes($hit->{ 'gtf.attributes' });
                    my $tid = exists $atts->{transcript_id} ? $atts->{transcript_id}[0] : undef;
                    if ($src eq 'refseq_canon' && $tid && exists($kgxref->{ $tid })) {
                        $features{description}{$kgxref->{ $tid }{description}}++;
                    }
                    $features{features}{$src}{ $atts->{gene_id}[0] }++;
                }
            }
        }
    }
    
    # if there are any features, they map to something (we consider these 'genic')
    if (exists $features{features}) {
        $str = join("\t",
            $buffer->{'bam.qname'}, # READ
            $seq_data->{sequence},   # SEQ
            length($seq_data->{sequence}), # LENGTH
            $seq_data->{hits},   # NUM_ALNS
            join(';', sort keys %{$features{location}}),
            (map { exists $features{features}{$_} ? join(';', keys %{$features{features}{$_}}) : '' } @VALID_SOURCES),
            exists $features{description} ? join(';', keys %{$features{description}} ) : '', 
             );
    }
    # no features, we consider these 'nongenic'
    else {  
        $str = join("\t",
        $buffer->{'bam.qname'}, # READ
        $seq_data->{sequence},   # SEQ
        length($seq_data->{sequence}), # LENGTH
        $seq_data->{hits},   # NUM_ALNS
        join(';', sort keys %{$features{location}})
         );
    }
    
    ($str, scalar keys %{$features{features}});
}
