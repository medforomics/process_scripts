#!/bin/bash
#germline_vc.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Path to Reference Genome with the file genome.fa"
  echo "-p  --Prefix for output file name"
  echo "-a  --Algorithm/Command: gatk, mpileup, speedseq, platypus "
  echo "Example: bash hisat.sh -p prefix -r /path/GRCh38 -a gatk"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:a:b:p:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
        a) algo=$OPTARG;;
        h) usage;;
    esac
done
function join_by { local IFS="$1"; shift; echo "$*"; }

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $index_path ]]; then
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
    dict="${index_path}/genome.dict"
else 
    echo "Missing Fasta File: ${index_path}/genome.fa"
    usage
fi
module load python/2.7.x-anaconda samtools/1.6 bedtools/2.26.0 snpeff/4.3q vcftools/0.1.14 parallel

for i in *.bam; do
    if [[ ! -f ${i}.bai ]]
    then
	samtools index --threads $SLURM_CPUS_ON_NODE $i
    fi
done

if [[ $algo == 'mpileup' ]]
then
    cut -f 1 ${index_path}/genomefile.5M.txt | parallel --delay 2 -j $SLURM_CPUS_ON_NODE "samtools mpileup -t 'AD,DP,INFO/AD' -ug -Q20 -C50 -r {} -f ${reffa} *.bam | bcftools call -vmO z -o ${pair_id}.samchr.{}.vcf.gz"
    vcf-concat ${pair_id}.samchr.*.vcf.gz |vcf-annotate -n --fill-type | bcftools norm -c s -f ${reffa} -w 10 -O v -o sam.vcf -
    java -jar $PICARD/picard.jar SortVcf I=sam.vcf O=${pair_id}.sam.vcf R=${reffa} CREATE_INDEX=TRUE
    bgzip ${pair_id}.sam.vcf
elif [[ $algo == 'hotspot' ]]
then
    samtools mpileup -d 99999 -t 'AD,DP,INFO/AD' -uf ${reffa} *.bam > ${pair_id}.mpi
    bcftools filter -i "AD[1]/DP > 0.01" ${pair_id}.mpi | bcftools filter -i "DP > 50" | bcftools call -m -A |vcf-annotate -n --fill-type |  bcftools norm -c s -f ${reffa} -w 10 -O z -o ${pair_id}.lowfreq.vcf.gz -
    java -jar $SNPEFF_HOME/SnpSift.jar annotate ${index_path}/cosmic.vcf.gz ${pair_id}.lowfreq.vcf.gz | java -jar $SNPEFF_HOME/SnpSift.jar filter "(CNT[*] >0)" - |bgzip > ${pair_id}.hotspot.vcf.gz
elif [[ $algo == 'speedseq' ]]
then
    module load speedseq/20160506
    speedseq var -t $SLURM_CPUS_ON_NODE -o ssvar ${reffa} *.bam
    vcf-annotate -n --fill-type ssvar.vcf.gz| bcftools norm -c s -f ${reffa} -w 10 -O z -o ${pair_id}.ssvar.vcf.gz -
elif [[ $algo == 'gatk' ]]
then
    module load gatk/3.7
    gvcflist=''
    for i in *.bam; do
	cut -f 1 ${index_path}/genomefile.5M.txt | parallel --delay 2 -j 10 "java -Djava.io.tmpdir=./ -Xmx32g -jar $GATK_JAR -R ${reffa} -D ${dbsnp} -T HaplotypeCaller -stand_call_conf 10 -A FisherStrand -A QualByDepth -A VariantType -A DepthPerAlleleBySample -A HaplotypeScore -A AlleleBalance -variant_index_type LINEAR -variant_index_parameter 128000 --emitRefConfidence GVCF -I $i -o ${i}.{}.chr.gatk.g.vcf -nct 2 -L {}"
	vcf-concat ${i}.*.chr.gatk.g.vcf |vcf-sort > gatk.g.vcf
	rm ${i}.*.chr.gatk.g.vcf
	java -jar $PICARD/picard.jar SortVcf I=gatk.g.vcf O=${i}.gatk.g.vcf R=${reffa} CREATE_INDEX=TRUE
	gvcflist="$gvcflist --variant ${i}.gatk.g.vcf"
    done
    java -Djava.io.tmpdir=./ -Xmx32g -jar $GATK_JAR -R ${reffa} -D ${dbsnp} -T GenotypeGVCFs -o gatk.vcf $gvcflist
    bcftools norm -c s -f ${reffa} -w 10 -O v gatk.vcf | vcf-annotate -n --fill-type gatk.vcf | bgzip > ${pair_id}.gatk.vcf.gz
    tabix ${pair_id}.gatk.vcf.gz
elif [[ $algo == 'platypus' ]]
then
    module load platypus/gcc/0.8.1
    bamlist=`join_by , *.bam`
    Platypus.py callVariants --minMapQual=10 --mergeClusteredVariants=1 --nCPU=$SLURM_CPUS_ON_NODE --bamFiles=${bamlist} --refFile=${reffa} --output=platypus.vcf
    vcf-sort platypus.vcf |vcf-annotate -n --fill-type -n |bgzip > platypus.vcf.gz
    tabix platypus.vcf.gz
    bcftools norm -c s -f ${reffa} -w 10 -O z -o ${pair_id}.platypus.vcf.gz platypus.vcf.gz
fi
