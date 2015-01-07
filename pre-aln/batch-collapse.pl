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
    dir         => 'nfiltered_seqs',
    ext         => '.trimmed.nfiltered.fastq.gz',
    threads     => 4,
    jobs        => 6,
    name        => 'fastx_derep',
    queue       => 'biotech',
    #db          => undef,
    #scripts     => '/home/groups/hpcbio/dev/pyrimidine/hpcbio/smalheiser_scripts/v2/'
);

GetOptions( \%args,
    'dir:s',
    'ext:s',
    'threads:i',
    'jobs:i',
    'name:s',
    'queue:s',
    #'db:s',
    #'scripts:s'
);

my @files = <$args{dir}/*$args{ext}>;
my $fc = @files;

die "No files found" if !@files;

my $tmp_dir = 'tmp';
mkdir $tmp_dir if !-e $tmp_dir;

my $ts = localtime->mdy('').'_'.localtime->hour().localtime->min();

my $script_dir = File::Spec->catdir('.', $tmp_dir, "$args{name}_$ts");
mkdir $script_dir if !-e $script_dir;

my $file_list = File::Spec->catfile('.', $script_dir, 'trimmed_file_list');

open(my $listfh, '>', $file_list);
print $listfh join("\n", @files);
close $listfh;

my $shell_script = <<SCRIPT; # heredoc
#!/bin/bash
#PBS -S /bin/bash
#PBS -N $args{name}
#PBS -t 1-$fc%$args{jobs}
#PBS -l nodes=1:ppn=$args{threads}
#PBS -e $script_dir/$args{name}_$ts.e
#PBS -o $script_dir/$args{name}_$ts.o
#PBS -l mem=16gb
#PBS -q $args{queue}

module load seqtk
module load fastx_toolkit/0.0.13

cd \$PBS_O_WORKDIR

FASTQ=`cat $file_list | tail -n +\${PBS_ARRAYID} | head -1`

SAMPLE_NAME=\$( basename \$FASTQ $args{ext} )

if [ ! -d collapsed_seqs ]; then
    mkdir collapsed_seqs
fi

if [ ! -e "collapsed_seqs/\${SAMPLE_NAME}.derep.fasta" ]; then

    zcat \$FASTQ | fastx_collapser -Q33 -v -o ./collapsed_seqs/\$SAMPLE_NAME.derep.fasta 2> ./collapsed_seqs/\$SAMPLE_NAME_collapsed.log

fi

SCRIPT

my $tmp = File::Spec->catfile($script_dir, "$args{name}_$ts.sh");

open(my $scrfh, '>', $tmp);
print $scrfh $shell_script;
close $tmp;

my $jobid = `qsub $tmp`;

$jobid =~ s/(\d+\[\])\..*$/$1/;

#my $old_data = '';

#if (exists $args{db}) {
#    $old_data = File::Spec->catfile($script_dir, "old_sequences.fasta");
#    my $generate_old = <<SCRIPT; # heredoc
##!/bin/bash
##PBS -S /bin/bash
##PBS -N $args{name}
##PBS -l nodes=1:ppn=$args{threads}
##PBS -e $script_dir/$args{name}_${ts}_final.e
##PBS -o $script_dir/$args{name}_${ts}_final.o
##PBS -W depend=afterokarray:$jobid
##PBS -l mem=16gb
##PBS -q $args{queue}
#
#module load perl
#cd \$PBS_O_WORKDIR
#
#perl $args{scripts}/utils/dump_old_fasta.pl -d $args{db} > $old_data
#
#SCRIPT
#    $tmp = File::Spec->catfile($script_dir, "$args{name}_${ts}_final.sh");
#
#    open($scrfh, '>', $tmp);
#    print $scrfh $generate_old;
#    close $tmp;
#    
#    $jobid = `qsub $tmp`;
#}

my $collapse_all = <<SCRIPT; # heredoc
#!/bin/bash
#PBS -S /bin/bash
#PBS -N $args{name}
#PBS -l nodes=1:ppn=$args{threads}
#PBS -e $script_dir/$args{name}_${ts}_final.e
#PBS -o $script_dir/$args{name}_${ts}_final.o
#PBS -W depend=afterokarray:$jobid
#PBS -l mem=16gb
#PBS -q $args{queue}

module load seqtk
module load fastx_toolkit/0.0.13

cd \$PBS_O_WORKDIR

cat ./collapsed_seqs/*.derep.fasta | fastx_collapser -Q33 -v 2> ./collapsed_seqs/all_collapsed.log | perl -p -e 's/>(\\d+)-(\\d+)/>srna\$1 \$2/' > collapsed_seqs/all.collapsed.fasta

SCRIPT

$tmp = File::Spec->catfile($script_dir, "$args{name}_${ts}_final.sh");

open($scrfh, '>', $tmp);
print $scrfh $collapse_all;
close $tmp;

$jobid = `qsub $tmp`;

# we can to reconcile the 
#if ($args{db}) {
#    my $reconcile_old = <<SCRIPT; # heredoc
##!/bin/bash
##PBS -S /bin/bash
##PBS -N $args{name}
##PBS -l nodes=1:ppn=$args{threads}
##PBS -e $script_dir/$args{name}_${ts}_final.e
##PBS -o $script_dir/$args{name}_${ts}_final.o
##PBS -W depend=afterokarray:$jobid
##PBS -l mem=16gb
##PBS -q $args{queue}
#
#module load perl
#
#cd \$PBS_O_WORKDIR
#
#perl $args{scripts}/utils/reconcile_old.pl -d $args{db} > $old_data
#
#SCRIPT
#}

