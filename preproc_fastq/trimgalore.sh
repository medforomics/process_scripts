#!/bin/bash
#trimgalore.sh

pair_id=$1
read1=$2
read2=$3
r1base="${r1%.fastq*}"
r2base="${r2%.fastq*}"

module load trimgalore/0.4.1 cutadapt/1.9.1
if [ $R1 == $R2 ]
then
    trim_galore -q 25 --illumina --gzip --length 35 ${read1}
    mv ${r1base}_trimmed.fq.gz ${pair_id}.trim.R1.fastq.gz
    cp ${pair_id}.trim.R1.fastq.gz ${pair_id}.trim.R2.fastq.gz 
else
    trim_galore --paired -q 25 --illumina --gzip --length 35 ${read1} ${read2}
    mv ${r1base}_val_1.fq.gz ${pair_id}.trim.R1.fastq.gz
    mv ${r2base}_val_2.fq.gz ${pair_id}.trim.R2.fastq.gz
fi
