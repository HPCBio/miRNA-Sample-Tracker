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
    ext         => 'fastq.gz',
    threads     => 8,
    adaptors    => '/home/apps/trimmomatic/trimmomatic-0.30/adapters/TruSeq3-SE.fa',
    # TODO: allow for optional steps (e.g. an array)
    leading     => 30,
    trailing    => 30,
    minlen      => 15,
    seed_mm     => 2,
    pal_clip    => 30,
    sim_clip    => 10,
    name        => 'trimmomatic_trim',
    queue       => 'biotech'
);

GetOptions( \%args,
    'file:s',
    'ext:s',
    'dir:s',
    'threads:i',
    'adaptors:s',
    'leading:i',
    'trailing:i',
    'minlen:i',
    'seed_mm:i',
    'pal_clip:i',
    'sim_clip:s',
    'name:s',
    'queue:s'
);

my $bam_dir = shift;
my @files = <$args{dir}/*.$args{ext}>;
my $fc = @files;

my $tmp_dir = 'tmp';
mkdir $tmp_dir if !-e $tmp_dir;

my $file_list = File::Spec->catfile('.', $tmp_dir, 'file_list');
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
#PBS -t 1-$fc%10
#PBS -l nodes=1:ppn=$args{threads}
#PBS -e $scriptdir/${ts}_$args{name}.e
#PBS -o $scriptdir/${ts}_$args{name}.o
#PBS -q $args{queue}

module load trimmomatic/0.30
CLASSPATH=/home/apps/trimmomatic/trimmomatic-0.30

cd \$PBS_O_WORKDIR

FQ_R1=`cat $file_list | tail -n +\${PBS_ARRAYID} | head -1`

SAMPLE_R1=\$( basename \$FQ_R1 .$args{ext} )

if [ ! -d trimmed_seqs_test ]; then
    mkdir trimmed_seqs_test
fi

set -x

if [ ! -e "trimmed_seqs_test/\${SAMPLE_R1}.trimmed.fastq.gz" ]; then
    for STEP in \$( seq 2 3 ); do
        java -Xmx6g -jar \$CLASSPATH/trimmomatic-0.30.jar SE \\
            -threads \$PBS_NUM_PPN \\
            -phred33 \\
            -trimlog /state/partition1/\${SAMPLE_R1}.\$STEP.trimlog.txt \\
            \$FQ_R1 \\
            /state/partition1/\${SAMPLE_R1}.\$STEP.trimmed.fastq.gz \\
            ILLUMINACLIP:$args{adaptors}:$args{seed_mm}:$args{pal_clip}:\$STEP \\
            LEADING:$args{leading} \\
            TRAILING:$args{trailing} \\
            MINLEN:$args{minlen} 2> /state/partition1/\${SAMPLE_R1}.\$STEP.trimmomatic.summary.log

        cut -d ' ' -f 3 /state/partition1/\${SAMPLE_R1}.\$STEP.trimlog.txt | sort | uniq -c > /state/partition1/\${SAMPLE_R1}.\$STEP.trimmomatic.hist.txt

        pigz -9 -p \$PBS_NUM_PPN /state/partition1/\${SAMPLE_R1}.\$STEP.trimlog.txt

        mv /state/partition1/\${SAMPLE_R1}.\$STEP.* ./trimmed_seqs_test
    done
fi

SCRIPT

my $tmp = File::Spec->catfile($scriptdir, "${ts}_$args{name}.sh");

open(my $scrfh, '>', $tmp);
print $scrfh $shell_script;
close $tmp;

my $status = `qsub $tmp`;
print $status;
