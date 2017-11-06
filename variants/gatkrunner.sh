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

if [[ -z $SLURM_CPUS_ON_NODE ]]
then
    SLURM_CPUS_ON_NODE=1
fi
if [[ -a "${index_path}/dbSnp.vcf.gz" ]]
then
    dbsnp="${index_path}/dbSnp.vcf.gz"
else 
    echo "Missing dbSNP File: ${index_path}/dbSnp.vcf.gz"
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

module load gatk/3.7 samtools/1.6
samtools index -@ $SLURM_CPUS_ON_NODE ${sbam}

if [[ $algo == 'gatkbam_rna' ]]
then
    module load picard/2.10.3
    java -Xmx4g -jar $PICARD/picard.jar CleanSam INPUT=${sbam} OUTPUT=${pair_id}.clean.bam
    java -Xmx4g -jar $PICARD/picard.jar ReorderSam I=${pair_id}.clean.bam O=${pair_id}.sort.bam R=${reffa} CREATE_INDEX=TRUE 
    #java -Xmx4g -jar $PICARD/picard.jar AddOrReplaceReadGroups INPUT=${pair_id}.clean.bam O=${pair_id}.rg_added_sorted.bam SO=coordinate RGID=${pair_id} RGLB=tx RGPL=illumina RGPU=barcode RGSM=${pair_id}
    #samtools index -@ $SLURM_CPUS_ON_NODE ${pair_id}.clean.bam
    java -Xmx4g -jar $GATK_JAR -T SplitNCigarReads -R ${reffa} -I ${pair_id}.clean.bam -o ${pair_id}.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
    java -Xmx32g -jar $GATK_JAR -T RealignerTargetCreator -known ${knownindel} -R ${reffa} -o ${pair_id}.bam.list -I ${pair_id}.split.bam -nt $SLURM_CPUS_ON_NODE -nct 1
    java -Xmx32g -jar $GATK_JAR -I ${pair_id}.split.bam -R ${reffa} --filter_mismatching_base_and_quals -T IndelRealigner -targetIntervals ${pair_id}.bam.list -o ${pair_id}.realigned.bam
    java -Xmx32g -jar $GATK_JAR -l INFO -R ${reffa} --knownSites ${dbsnp} -I ${pair_id}.realigned.bam -T BaseRecalibrator -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o ${pair_id}.recal_data.grp -nt 1 -nct $SLURM_CPUS_ON_NODE
    java -Xmx32g -jar $GATK_JAR -T PrintReads -R ${reffa} -I ${pair_id}.realigned.bam -BQSR ${pair_id}.recal_data.grp -o ${pair_id}.final.bam -nt 1 -nct $SLURM_CPUS_ON_NODE
elif [[ $algo == 'gatkbam' ]]
then
  java -Xmx32g -jar $GATK_JAR -T RealignerTargetCreator -known ${knownindel} -R ${reffa} -o ${pair_id}.bam.list -I ${sbam} -nt $SLURM_CPUS_ON_NODE -nct 1
  java -Xmx32g -jar $GATK_JAR -I ${sbam} -R ${reffa} --filter_mismatching_base_and_quals -T IndelRealigner -targetIntervals ${pair_id}.bam.list -o ${pair_id}.realigned.bam
  java -Xmx32g -jar $GATK_JAR -l INFO -R ${reffa} --knownSites ${dbsnp} -I ${pair_id}.realigned.bam -T BaseRecalibrator -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o ${pair_id}.recal_data.grp -nt 1 -nct $SLURM_CPUS_ON_NODE
  java -Xmx32g -jar $GATK_JAR -T PrintReads -R ${reffa} -I ${pair_id}.realigned.bam -BQSR ${pair_id}.recal_data.grp -o ${pair_id}.final.bam -nt 1 -nct $SLURM_CPUS_ON_NODE
fi

