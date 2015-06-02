#!/usr/bin/perl -w
use 5.010;
use strict;
use warnings;
use Time::Piece;
use File::Temp;
use File::Spec;
use Getopt::Long;

# TODO: this needs to be genericized (maybe a templating system?) and allow for
# job arrays

my %args = (
    dir         => '.',
    ext         => 'trimmed.fastq.gz',
    threads     => 6,
    name        => 'nfilter_array',
    queue       => 'biotech'
);

GetOptions( \%args,
    'dir:s',
    'ext:s',
    'threads:i',
    'name:s',
    'queue:s'
);

#my $bam_dir = shift;
my @files = <$args{dir}/*.$args{ext}>;
my $fc = @files;

my $tmp_dir = 'tmp';
mkdir $tmp_dir if !-e $tmp_dir;

my $file_list = File::Spec->catfile('.', $tmp_dir, 'trimmed_file_list');
my $ts = localtime->mdy('').'_'.localtime->hour().localtime->min();

my $scriptdir = File::Spec->catdir('.', $tmp_dir, "${ts}_$args{name}");
mkdir $scriptdir if !-e $scriptdir;

open(my $listfh, '>', $file_list);
print $listfh join("\n", @files);
close $listfh;

my $shell_script = <<SCRIPT; # heredoc
#!/bin/bash
#PBS -S /bin/bash
#PBS -N $args{name}
#PBS -t 1-$fc%20
#PBS -l nodes=1:ppn=$args{threads}
#PBS -e $scriptdir/$args{name}.e
#PBS -o $scriptdir/$args{name}.o
#PBS -q $args{queue}

cd \$PBS_O_WORKDIR

set -x

FQ_R1=`cat $file_list | tail -n +\${PBS_ARRAYID} | head -1`

SAMPLE_R1=\$( basename \$FQ_R1 .$args{ext} )

if [ ! -d nfiltered_seqs ]; then
    mkdir nfiltered_seqs
fi

if [ ! -e "processed_seqs/\${SAMPLE_R1}.trimmed.nfiltered.fastq.gz" ]; then
    filter_ns.pl \$FQ_R1 \$SAMPLE_R1 > /state/partition1/\${SAMPLE_R1}.trimmed.nfiltered.fastq

    pigz -p \$PBS_NP /state/partition1/\${SAMPLE_R1}.trimmed.nfiltered.fastq

    mv /state/partition1/\${SAMPLE_R1}.* ./nfiltered_seqs
    mv \${SAMPLE_R1}.* ./nfiltered_seqs
fi

SCRIPT

my $tmp = File::Spec->catfile($scriptdir, "$args{name}.sh");

open(my $scrfh, '>', $tmp);
print $scrfh $shell_script;
close $tmp;

my $status = `qsub $tmp`;
print $status;
