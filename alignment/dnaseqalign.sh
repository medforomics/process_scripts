#!/bin/bash
#dnaseqalign.sh

usage() {
  echo "-h Help documentation for dnaseqalign.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-x  --FastQ R1"
  echo "-y  --FastQ R2"
  echo "-p  --Prefix for output file name"
  echo "-u UMI"
  echo "Example: bash dnaseqalign.sh -p prefix -u 1 -r /project/shared/bicf_workflow_ref/GRCh38 -x SRR1551047_1.fastq.gz  -y SRR1551047_2.fastq.gz"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:x:y:p:uh opt
do
    case $opt in
        r) index_path=$OPTARG;;
        x) fq1=$OPTARG;;
        y) fq2=$OPTARG;;
	u) umi='umi';;
        p) pair_id=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $fq1 ]]; then
    usage
fi

if [[ -z $SLURM_CPUS_ON_NODE ]]
then
    SLURM_CPUS_ON_NODE=1
fi

module load bwakit/0.7.15 bwa/intel/0.7.15 samtools/1.6 picard/2.10.3

baseDir="`dirname \"$0\"`"

if [[ $fq1 == $fq2 ]]
then
    bwa mem -M -t $SLURM_CPUS_ON_NODE -R "@RG\tID:${pair_id}\tLB:tx\tPL:illumina\tPU:barcode\tSM:${pair_id}" ${index_path}/genome.fa ${fq1} > out.sam
else
    bwa mem -M -t $SLURM_CPUS_ON_NODE -R "@RG\tID:${pair_id}\tLB:tx\tPL:illumina\tPU:barcode\tSM:${pair_id}" ${index_path}/genome.fa ${fq1} ${fq2} > out.sam
fi

if [[ $umi == 'umi' ]]
then
    k8 /cm/shared/apps/bwakit/0.7.15/bwa-postalt.js -p ${pair_id}.hla ${index_path}/genome.fa.alt out.sam | python ${baseDir}/add_umi_sam.py -s - -o output.unsort.bam
elif [[ $index_path == '/project/shared/bicf_workflow_ref/GRCh38' ]]
then
    k8 /cm/shared/apps/bwakit/0.7.15/bwa-postalt.js -p ${pair_id}.hla ${index_path}/genome.fa.alt out.sam| samtools view -1 - > output.unsort.bam
else 
    samtools view -1 -o output.unsort.bam out.sam
fi

samtools sort -n --threads $SLURM_CPUS_ON_NODE -o output.dups.bam output.unsort.bam
java -Djava.io.tmpdir=./ -Xmx4g  -jar $PICARD/picard.jar FixMateInformation ASSUME_SORTED=TRUE SORT_ORDER=coordinate ADD_MATE_CIGAR=TRUE I=output.dups.bam O=${pair_id}.bam
samtools index -@ $SLURM_CPUS_ON_NODE ${pair_id}.bam
