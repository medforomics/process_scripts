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

source /etc/profile.d/modules.sh
user=$USER
module load gatk/4.x singularity/2.6.1
mkdir /tmp/${user}
export TMP_HOME=/tmp/${user}

samtools index -@ $SLURM_CPUS_ON_NODE ${sbam}

if [[ $algo == 'gatkbam_rna' ]]
then
    module load picard/2.10.3
    java -Xmx4g -jar $PICARD/picard.jar CleanSam INPUT=${sbam} OUTPUT=${pair_id}.clean.bam
    java -Xmx4g -jar $PICARD/picard.jar ReorderSam I=${pair_id}.clean.bam O=${pair_id}.sort.bam R=${reffa} CREATE_INDEX=TRUE 
    java -Xmx4g -jar $PICARD/picard.jar AddOrReplaceReadGroups INPUT=${pair_id}.clean.bam O=${pair_id}.rg_added_sorted.bam SO=coordinate RGID=${pair_id} RGLB=tx RGPL=illumina RGPU=barcode RGSM=${pair_id}
    samtools index -@ $SLURM_CPUS_ON_NODE ${pair_id}.clean.bam
    java -Xmx4g -jar $GATK_JAR -L ${index_path}/../gatk_regions.list -T SplitNCigarReads -R ${reffa} -I ${pair_id}.sort.bam -o ${pair_id}.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
    java -Xmx32g -jar $GATK_JAR -L ${index_path}/../gatk_regions.list -T RealignerTargetCreator -known ${knownindel} -R ${reffa} -o ${pair_id}.bam.list -I ${pair_id}.split.bam -nt 8 -nct 1
    java -Xmx16g -jar $GATK_JAR -L ${index_path}/../gatk_regions.list -I ${pair_id}.split.bam -R ${reffa} --filter_mismatching_base_and_quals -T IndelRealigner -targetIntervals ${pair_id}.bam.list -o ${pair_id}.realigned.bam
    java -Xmx16g -jar $GATK_JAR -l INFO -R ${reffa} --knownSites ${dbsnp} -I ${pair_id}.realigned.bam -T BaseRecalibrator -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o ${pair_id}.recal_data.grp -nt 1 -nct 8
    java -Xmx16g -jar $GATK_JAR -L ${index_path}/../gatk_regions.list -T PrintReads -R ${reffa} -I ${pair_id}.realigned.bam -BQSR ${pair_id}.recal_data.grp -o ${pair_id}.final.bam -nt 1 -nct 8

elif [[ $algo == 'gatkbam' ]]
then
    singularity exec -H /tmp/${user} /project/apps/singularity-images/gatk4/gatk-4.x.simg /gatk/gatk --java-options "-Xmx32g" BaseRecalibrator -I ${i} --known-sites ${gatk4_dbsnp} -R ${reffa} -O ${prefix}.recal_data.table --use-original-qualities
    singularity exec -H /tmp/${user} /project/apps/singularity-images/gatk4/gatk-4.x.simg /gatk/gatk --java-options "-Xmx32g" ApplyBQSR -I ${i} -R ${reffa} -O ${prefix}.final.bam --use-original-qualities -bqsr ${prefix}.recal_data.table
    samtools index -@ $SLURM_CPUS_ON_NODE ${pair_id}.final.bam

elif [[ $algo == 'abra2' ]]
then
  module load abra2/2.18
  mkdir tmpdir
  java  -Xmx16G -jar /cm/shared/apps/abra2/lib/abra2.jar --in ${sbam}  --in-vcf /archive/PHG/PHG_Clinical/phg_workflow/analysis/awesomeproject/GoldIndels.vcf --out ${pair_id}.final.bam --ref ${reffa} --threads $SLURM_CPUS_ON_NODE --tmpdir tmpdir
  samtools index -@ $SLURM_CPUS_ON_NODE ${pair_id}.final.bam
fi
