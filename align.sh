#!/bin/bash
#PBS -S /bin/bash
#PBS -l mem=24gb
#PBS -A Smalheiser
#PBS -q biotech
#PBS -l nodes=1:ppn=8

####################################################################
#
# Must pass in via 'qsub -v FQ_R1=r1.fq,FQ_R2=r2.fq,INDEX=foo.ndx align.sh'
#
# INDEX = Novoalign index file
# SEQUENCE = R1 reads to be aligned (these are expected to be FASTA)
#
####################################################################

# add Novoalign and samtools to PATH

module load novocraft/3.02
module load samtools

cd $PBS_O_WORKDIR

SAMPLE=$( basename $( readlink -e $SEQUENCE ) )
SAMPLE=${SAMPLE%%.*}
DB_GPFS=$( basename ${INDEX%%.*} )

DB_GPFS=${DB_GPFS}_ALL

DB_BASE=/state/partition1/${DB_GPFS}_ALL

if [ ! -d $DB_BASE ]; then
    mkdir -p $DB_BASE
fi

if [ ! -d $DB_GPFS ]; then
    mkdir -p $DB_GPFS
fi

# -c 4 -d genome_index/ath_miRBase_19.nix
# -f Step3_trimmed_fastq_file/1_ATCACG_L001_R1_001_trimed_longer_15.fastq
# -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC -l 15 -r ALL >Step5_alignment_file/1_ATCACG_ath_new_alignment.txt

if [ ! -f $DB_GPFS/$SAMPLE.bam ]; then

    novoalign -c $( expr $PBS_NUM_PPN - 1 )  \
        -d $INDEX \
        -f $SEQUENCE \
        -l 15 \
        -r "${MODE:=ALL}" \
        -o SAM \
        -m \
        2> $DB_BASE/$SAMPLE.run.log | samtools view -bS - > $DB_BASE/$SAMPLE.bam

    if [ ! -f $SAMPLE.sorted.bam ]; then
        # Sort BAM
        novosort -c $PBS_NUM_PPN -t /state/partition1 -m 20G -i $DB_BASE/$SAMPLE.bam -o $DB_BASE/$SAMPLE.sorted.bam
        samtools flagstat $DB_BASE/$SAMPLE.sorted.bam > $DB_BASE/$SAMPLE.flagstats.txt
        samtools sort -n -@ $PBS_NUM_PPN $DB_BASE/$SAMPLE.bam $DB_BASE/$SAMPLE.name
        novosort -n -c $PBS_NUM_PPN -t /state/partition1 -m 20G -i $DB_BASE/$SAMPLE.bam -o $DB_BASE/$SAMPLE.name.bam
    fi

    mv $DB_BASE/${SAMPLE}* $DB_GPFS/

    # this will fail if there is an error
    rmdir $DB_BASE/

fi

