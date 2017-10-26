#!/bin/bash
#hisat.sh

usage() {
  echo "-h Help documentation for hisat.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-a  --FastQ R1"
  echo "-b  --FastQ R2"
  echo "-p  --Prefix for output file name"
  echo "Example: bash hisat.sh -p prefix -r GRCh38 -a SRR1551047_1.fastq.gz  -b SRR1551047_2.fastq.gz"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:a:b:p:h opt
do
    case $opt in
        r) refgeno=$OPTARG;;
        a) fq1=$OPTARG;;
        b) fq2=$OPTARG;;
        p) pair_id=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $fq1 ]]; then
    usage
fi

if [ $refgeno == 'GRCh38' ] || [ $refgeno == 'GRCm38' ]; then
    index_path=/project/shared/bicf_workflow_ref/${refgeno}/hisat_index/
else
    usage
fi

if [[ -z $SLURM_CPUS_ON_NODE ]]
then
    SLURM_CPUS_ON_NODE=1
fi
module load hisat2/2.1.0-intel samtools/gcc/1.6 picard/2.10.3
if [ $fq1 == $fq2 ]
then
    hisat2 -p $SLURM_CPUS_ON_NODE --rg-id ${pair_id} --rg LB:tx --rg PL:illumina --rg PU:barcode --rg SM:${pair_id} --add-chrname --no-unal --dta -x ${index_path}/genome -U ${fq1} -S out.sam --summary-file ${pair_id}.alignerout.txt
else
    hisat2 -p $SLURM_CPUS_ON_NODE --rg-id ${pair_id} --rg LB:tx --rg PL:illumina --rg PU:barcode --rg SM:${pair_id} --add-chrname --no-unal --dta -x ${index_path}/genome -1 ${fq1} -2 ${fq2} -S out.sam --summary-file ${pair_id}.alignerout.txt
fi
samtools view -1 --threads $SLURM_CPUS_ON_NODE -o output.bam out.sam
samtools sort -@ $SLURM_CPUS_ON_NODE -O BAM -n -o  output.nsort.bam output.bam
java -jar $PICARD/picard.jar FixMateInformation ASSUME_SORTED=TRUE SORT_ORDER=coordinate ADD_MATE_CIGAR=TRUE I=output.nsort.bam O=${pair_id}.bam
samtools index -@ $SLURM_CPUS_ON_NODE ${pair_id}.bam
