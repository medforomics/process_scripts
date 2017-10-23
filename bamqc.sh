#!/bin/bash
#trimgalore.sh

pair_id=$1
bam=$2

module load samtools/1.4.1 fastqc/0.11.5
samtools flagstat ${bam} > ${pair_id}.flagstat.txt
fastqc -f bam ${bam}
