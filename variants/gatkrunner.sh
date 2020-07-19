#!/bin/bash
#gatkrunner.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-b  --BAM File"
  echo "-p  --Prefix for output file name"
  echo "-a  --Algorithm/Command"
  echo "Example: bash hisat.sh -p prefix -r GRCh38 -b File.bam"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:a:b:p:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        b) sbam=$OPTARG;;
        p) pair_id=$OPTARG;;
        a) algo=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $sbam ]] || [[ -z $index_path ]]
then
    usage
fi

NPROC=$SLURM_CPUS_ON_NODE
if [[ -z $NPROC ]]
then
    NPROC=`nproc`
fi
if [[ -a "${index_path}/dbSnp.gatk4.vcf.gz" ]]
then
    dbsnp="${index_path}/dbSnp.gatk4.vcf.gz"
else 
    echo "Missing dbSNP File: ${index_path}/dbSnp.gatk4.vcf.gz"
    usage
fi
if [[ -a "${index_path}/GoldIndels.vcf.gz" ]]
then
    knownindel="${index_path}/GoldIndels.vcf.gz"
else 
    echo "Missing InDel File: ${index_path}/GoldIndels.vcf.gz"
    usage
fi
if [[ -a "${index_path}/genome.fa" ]]
then
    reffa="${index_path}/genome.fa"
else 
    echo "Missing Fasta File: ${index_path}/genome.fa"
    usage
fi

source /etc/profile.d/modules.sh
module load gatk/4.1.4.0 samtools/gcc/1.8
which samtools
samtools index -@ $NPROC ${sbam}

gatk --java-options "-Xmx32g" BaseRecalibrator -I ${sbam} --known-sites ${index_path}/dbSnp.gatk4.vcf.gz -R ${reffa} -O ${pair_id}.recal_data.table --use-original-qualities
gatk --java-options "-Xmx32g" ApplyBQSR -I ${sbam} -R ${reffa} -O ${pair_id}.final.bam --use-original-qualities -bqsr ${pair_id}.recal_data.table
samtools index -@ $NPROC ${pair_id}.final.bam
