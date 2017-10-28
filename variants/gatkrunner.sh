#!/bin/bash
#gatkrunner.sh

usage() {
  echo "-h Help documentation for hisat.sh"
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
        r) refgeno=$OPTARG;;
        b) sbam=$OPTARG;;
        p) pair_id=$OPTARG;;
        a) algo=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $bam ] || [[ -z $refgeno ]]]; then
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

dbsnp="${index_path}/Snp.vcf.gz"
knownindel="${index_path}/Snp.vcf.gz"
reffa="${index_path}/genome.fa"

module load gatk/3.7 samtools/1.6

if [[ $algo == 'gatkbam' ]]
then
    samtools index ${sbam}
  java -Xmx32g -jar $GATK_JAR -T RealignerTargetCreator -known ${knownindel} -R ${reffa} -o ${pair_id}.bam.list -I ${sbam} -nt $SLURM_CPUS_ON_NODE -nct 1
  java -Xmx32g -jar $GATK_JAR -I ${sbam} -R ${reffa} --filter_mismatching_base_and_quals -T IndelRealigner -targetIntervals ${pair_id}.bam.list -o ${pair_id}.realigned.bam
  java -Xmx32g -jar $GATK_JAR -l INFO -R ${reffa} --knownSites ${dbsnp} -I ${pair_id}.realigned.bam -T BaseRecalibrator -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o ${pair_id}.recal_data.grp -nt 1 -nct $SLURM_CPUS_ON_NODE
  java -Xmx32g -jar $GATK_JAR -T PrintReads -R ${reffa} -I ${pair_id}.realigned.bam -BQSR ${pair_id}.recal_data.grp -o ${pair_id}.final.bam -nt 1 -nct 8
elif
then

else

fi
