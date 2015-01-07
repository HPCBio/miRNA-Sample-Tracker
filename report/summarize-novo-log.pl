#!/usr/bin/env perl
use 5.016;
use strict;
use warnings;
use autodie;
use Getopt::Long;
use File::Basename;
use File::Spec;
use Data::Dumper;
use File::Find;

my %data;

my %args = (
    dir         => '.',
    ext         => '.run.log',
);

GetOptions( \%args,
    'dir:s',
    'ext:s',
);

my %keys;

my %primary = (
    'Read Sequences'        => 1,
    'Aligned'               => 2,
    'Unique Alignment'      => 3,
    'Gapped Alignment'      => 4,
    'Homopolymer Filter'    => 5,
    'Quality Filter'        => 6,
);

my @files; # = map { File::Spec->rel2abs($_) } <$args{dir}/*$args{ext}>;

my $cb = sub {
    if ($_ =~ /$args{ext}/) {
        push @files, $File::Find::name;
    }
};

find($cb, $args{dir});

for my $file (@files) {
    #(my $sample = File::Basename::basename($file)) =~ s/_[ATGC]{6}_L\d+_(R[12]).*/$1/;
    open(my $LOG, '<', $file);

    while (<$LOG>) {
        chomp;
        s/^#\s+//;
        my ($k, $v);
        if (/^novoalign/) {
            if (/(V\d+\.\d+\.\d+)/) {
                ($k, $v) = ('version', $1);
            } else {
                ($k, $v) = ('command-line',  $_);
            }
        } elsif (/(Starting|Done)\s+at\s+(.*)/) {
            ($k, $v) = ("$1 Time", $2);
        } elsif (/:/) {
            ($k, $v) = split(/\s*:\s*/, $_, 2);
        }
        next unless defined $k;
        $data{$file}{$k} = $v;
        $keys{$k}++;
    }

    close $LOG;
}

# warning: this is an ugly hack!
for my $k (sort keys %keys) {
    $primary{$k} = 999 if !exists $primary{$k};
}

my @order = sort {
    $primary{$a} <=> $primary{$b}
    } keys %primary;

say join("\t", "Sample", @order);

for my $file (sort keys %data) {
    #(my $nm = File::Basename::basename($file)) =~ s/\.trimmed\.run\.log//;
    say join("\t", $file, map {$data{$file}{$_}} @order);
}

__END__

# novoalign (V3.02.00 - Build Oct 22 2013 @ 08:35:33 - A short read aligner with qualities.
# (C) 2008,2009,2010,2011 NovoCraft Technologies Sdn Bhd.
# License file: /home/apps/novocraft/novocraft-3.02/novoalign.lic
# Licensed to University of Illinois
#  novoalign -c 7 -d /home/groups/hpcbio/projects/barnes-singleton-wani/April-2014/databases/A506.fasta.nix -f /home/groups/hpcbio/projects/barnes-singleton-wani/April-2014/processed_seqs/4C_I_CCGTCC_L007_R1_001.trimmed.fastq.gz -l 15 -r RANDOM -o SAM -# 1M -m
# Starting at Tue Apr 15 17:00:51 2014
# Interpreting input files as Illumina FASTQ, Casava Pipeline 1.8 & later.
# Index Build Version: 3.2
# Hash length: 10
# Step size: 1
#     Read Sequences:  1000000
#            Aligned:   418000
#   Unique Alignment:   414180
#   Gapped Alignment:    55493
#     Quality Filter:        8
# Homopolymer Filter:     7686
#       Elapsed Time: 142.692 (secs.)
#           CPU Time: 16.60 (min.)
# Done at Tue Apr 15 17:03:14 2014
