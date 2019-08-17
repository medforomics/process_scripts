#!/bin/bash
#svcalling.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Path to Reference Genome with the file genome.fa"
  echo "-p  --Prefix for output file name"
  echo "Example: bash svcalling.sh -p prefix -r /path/GRCh38 -a gatk"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:b:l:p:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        b) sbam=$OPTARG;;	
        p) pair_id=$OPTARG;;
	l) idtbed=$OPTARG;;
        h) usage;;
    esac
done
function join_by { local IFS="$1"; shift; echo "$*"; }

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"


# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $index_path ]]; then
    usage
fi
if [[ -z $SLURM_CPUS_ON_NODE ]]
then
    SLURM_CPUS_ON_NODE=1
fi

if [[ -a "${index_path}/genome.fa" ]]
then
    reffa="${index_path}/genome.fa"
    dict="${index_path}/genome.dict"
else 
    echo "Missing Fasta File: ${index_path}/genome.fa"
    usage

fi

source /etc/profile.d/modules.sh	


module load samtools/1.6 snpeff/4.3q vcftools/0.1.14 bcftools/gcc/1.8 bedtools/2.26.0 
stexe=`which samtools`

#flt3='chr13:28003274-28100592'

samtools view -@ $SLURM_CPUS_ON_NODE -L ${idtbed} ${sbam} | /project/shared/bicf_workflow_ref/seqprg/itdseek-1.2/itdseek.pl --refseq ${reffa} --samtools ${stexe} --bam ${sbam} | vcf-sort |bgzip > ${pair_id}.idtseek.vcf.gz
tabix ${pair_id}.idtseek.vcf.gz
bedtools intersect -header -b ${idtbed} -a ${pair_id}.idtseek.vcf.gz | java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config GRCh38.86 - |bgzip > ${pair_id}.idtseek_tandemdup.vcf.gz
