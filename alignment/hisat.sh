#!/bin/bash
#hisat.sh

pair_id=$1
index_path=$2
fq1=$3
fq2=$4
if [[ -z $SLURM_CPUS_ON_NODE ]]
then
    SLURM_CPUS_ON_NODE=1
fi

module load hisat2/2.1.0-intel samtools/1.4.1
if [ $fq1 == $fq2 ]
then
    hisat2 -p $SLURM_CPUS_ON_NODE --rg-id ${pair_id} --rg '@RG\\tID:${pair_id}\\tLB:tx\\tPL:illumina\\tPU:barcode\\tSM:${pair_id}' --add-chrname --no-unal --dta -x ${index_path}/genome -U ${fq1} -S out.sam --summary-file ${pair_id}.alignerout.txt
else
    hisat2 -p $SLURM_CPUS_ON_NODE --rg-id ${pair_id} --rg '@RG\\tID:${pair_id}\\tLB:tx\\tPL:illumina\\tPU:barcode\\tSM:${pair_id}' --add-chrname --no-unal --dta -x ${index_path}/genome -1 ${fq1} -2 ${fq2} -S out.sam --summary-file ${pair_id}.alignerout.txt
fi
samtools view -1 --threads $SLURM_CPUS_ON_NODE -o output.bam out.sam
#fixmateinfomation
samtools sort --threads $SLURM_CPUS_ON_NODE -o ${pair_id}.bam output.bam
