#!/bin/bash
#dnaseqalign.sh

usage() {
  echo "-h Help documentation for hisat.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-x  --FastQ R1"
  echo "-y  --FastQ R2"
  echo "-p  --Prefix for output file name"
  echo "Example: bash hisat.sh -p prefix -r GRCh38 -x SRR1551047_1.fastq.gz  -y SRR1551047_2.fastq.gz"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:a:x:y:p:h opt
do
    case $opt in
        r) refgeno=$OPTARG;;
        x) fq1=$OPTARG;;
        y) fq2=$OPTARG;;
	a) algo=$OPTARG;;
        p) pair_id=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $fq1 ]]; then
    usage
fi

if [ $refgeno == 'GRCh38' ] || [ $refgeno == 'GRCm38' ] || [ $refgeno == 'GRCh37' ]; then
    index_path=/project/shared/bicf_workflow_ref/${refgeno}
else
    usage
fi

if [[ -z $SLURM_CPUS_ON_NODE ]]
then
    SLURM_CPUS_ON_NODE=1
fi

module load bwakit/0.7.15 bwa/intel/0.7.15 samtools/1.6 bcftools/1.6 htslib/1.6 picard/2.10.3 speedseq/20160506
if ($fq1 != $fq2)
then
    bwa mem -M -t $SLURM_CPUS_ON_NODE -R "@RG\tID:${pair_id}\tLB:tx\tPL:illumina\tPU:barcode\tSM:${pair_id}" ${index_path}/genome.fa ${fq1} > out.sam
else
    bwa mem -M -t $SLURM_CPUS_ON_NODE -R "@RG\tID:${pair_id}\tLB:tx\tPL:illumina\tPU:barcode\tSM:${pair_id}" ${index_path}/genome.fa ${fq1} ${fq2} > out.sam
fi
if [ $refgeno == 'GRCh38' ]
then
    cat out.sam | k8 /cm/shared/apps/bwakit/0.7.15/bwa-postalt.js -p ${pair_id}.hla ${index_path}/genome.fa.alt | samtools view -1 - > output.unsort.bam
    run-HLA ${pair_id}.hla > ${pair_id}.hla.top 2> ${pair_id}.log.hla
    touch ${pair_id}.hla.HLA-dummy.gt
    cat ${pair_id}.hla.HLA*.gt | grep ^GT | cut -f2- > ${pair_id}.hla.all
else 
    samtools view -1 -o output.unsort.bam out.sam
fi 
samtools sort --threads $SLURM_CPUS_ON_NODE -o output.dups.bam output.unsort.bam
samtools view -h output.unsort.bam | samblaster -M -a --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 -d discordants.sam -s splitters.sam > temp.sam
gawk '{ if (\$0~"^@") { print; next } else { \$10="*"; \$11="*"; print } }' OFS="\\t" splitters.sam | samtools  view -S -b - | samtools sort -o ${pair_id}.splitters.bam -
gawk '{ if (\$0~"^@") { print; next } else { \$10="*"; \$11="*"; print } }' OFS="\\t" discordants.sam | samtools  view -S  -b - | samtools sort -o ${pair_id}.discordants.bam -
java -Djava.io.tmpdir=./ -Xmx4g  -jar \$PICARD/picard.jar FixMateInformation ASSUME_SORTED=TRUE SORT_ORDER=coordinate ADD_MATE_CIGAR=TRUE I=output.dups.bam O=${pair_id}.bam
samtools index -@ $SLURM_CPUS_ON_NODE ${pair_id}.bam
