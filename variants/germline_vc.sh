#!/bin/bash
#germline_vc.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Path to Reference Genome with the file genome.fa"
  echo "-p  --Prefix for output file name"
  echo "-a  --Algorithm/Command: gatk, mpileup, speedseq, platypus"
  echo "-t  --RNASeq Data"
  echo "Example: bash hisat.sh -p prefix -r /path/GRCh38 -a gatk"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:a:b:q:p:th opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
        a) algo=$OPTARG;;
	t) rna=1;;
	b) tbed=$OPTARG;;
	q) pon==$OPTARG;; 
        h) usage;;
    esac
done
function join_by { local IFS="$1"; shift; echo "$*"; }

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $index_path ]]; then
    usage
fi
NPROC=$SLURM_CPUS_ON_NODE
if [[ -z $NPROC ]]
then
    NPROC=`nproc`
fi
if [[ -s "${index_path}/dbSnp.vcf.gz" ]]
then
    dbsnp="${index_path}/dbSnp.vcf.gz"
    gatk4_dbsnp="${index_path}/dbSnp.gatk4.vcf.gz"
else 
    echo "Missing dbSNP File: ${index_path}/dbSnp.vcf.gz"
    usage
fi
if [[ -s "${index_path}/GoldIndels.vcf.gz" ]]
then
    knownindel="${index_path}/GoldIndels.vcf.gz"
else 
    echo "Missing InDel File: ${index_path}/GoldIndels.vcf.gz"
    usage
fi
if [[ -s "${index_path}/genome.fa" ]]
then
    reffa="${index_path}/genome.fa"
    dict="${index_path}/genome.dict"
else 
    echo "Missing Fasta File: ${index_path}/genome.fa"
    usage
fi
if [[ -f $pon ]]
then
    ponopt="--pon $pon"
else
    ponopt='';
fi

source /etc/profile.d/modules.sh
module load python/2.7.x-anaconda picard/2.10.3 samtools/gcc/1.8 bcftools/gcc/1.8 bedtools/2.26.0 snpeff/4.3q vcftools/0.1.14 parallel

for i in *.bam; do
    if [[ ! -f ${i}.bai ]]
    then
	samtools index -@ $NPROC $i
    fi
done

if [[ $algo == 'mpileup' ]]
then
    threads=`expr $NPROC - 10`
    bcftools mpileup --threads $threads -a 'INFO/AD,INFO/ADF,INFO/ADR,FORMAT/DP,FORMAT/SP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR' -Ou -A -d 1000000 -C50 -f ${reffa} *.bam | bcftools call -A --threads 10 -vmO z -o ${pair_id}.vcf.gz
    vcf-annotate -n --fill-type ${pair_id}.vcf.gz | bcftools norm -c s -f ${reffa} -w 10 -O v -o sam.vcf
    java -jar $PICARD/picard.jar SortVcf I=sam.vcf O=${pair_id}.sam.vcf R=${reffa} CREATE_INDEX=TRUE
    bgzip ${pair_id}.sam.vcf
elif [[ $algo == 'fb' ]]
then
    module load freebayes/gcc/1.2.0 parallel/20150122
    bamlist=''
    for i in *.bam; do
    bamlist="$bamlist --bam ${PWD}/${i}"
    done
    cut -f 1 ${index_path}/genomefile.5M.txt | parallel --delay 2 -j $NPROC "freebayes -f ${index_path}/genome.fa  --min-mapping-quality 0 --min-base-quality 20 --min-coverage 10 --min-alternate-fraction 0.01 -C 3 --use-best-n-alleles 3 -r {} ${bamlist} > fb.{}.vcf"
    vcf-concat fb.*.vcf | vcf-sort | vcf-annotate -n --fill-type | bcftools norm -c s -f ${reffa} -w 10 -O z -o ${pair_id}.fb.vcf.gz -
elif [[ $algo == 'platypus' ]]
then
    module load platypus/gcc/0.8.1
    bamlist=`join_by , *.bam`
    Platypus.py callVariants --minMapQual=0 --minReads=3 --mergeClusteredVariants=1 --nCPU=$NPROC --bamFiles=${bamlist} --refFile=${reffa} --output=platypus.vcf
    vcf-sort platypus.vcf |vcf-annotate -n --fill-type -n |bgzip > platypus.vcf.gz
    tabix platypus.vcf.gz
    bcftools norm -c s -f ${reffa} -w 10 -O z -o ${pair_id}.platypus.vcf.gz platypus.vcf.gz
elif [[ $algo == 'gatk' ]]
then
    user=$USER
    module load gatk/4.1.4.0
    gvcflist=''
    for i in *.bam; do
	prefix="${i%.bam}"
	echo ${prefix}
	gatk --java-options "-Xmx32g" HaplotypeCaller -R ${reffa} -I ${i} -A FisherStrand -A QualByDepth  -A DepthPerAlleleBySample -A TandemRepeat --emit-ref-confidence GVCF -O haplotypecaller.vcf.gz
	java -jar $PICARD/picard.jar SortVcf I=haplotypecaller.vcf.gz O=${prefix}.gatk.g.vcf R=${reffa} CREATE_INDEX=TRUE
	
	gvcflist="$gvcflist -V ${prefix}.gatk.g.vcf"
    done
    interval=`cat ${reffa}.fai |cut -f 1 |grep -v decoy |grep -v 'HLA' |grep -v alt |grep -v 'chrUn' |grep -v 'random' | perl -pe 's/\n/ -L /g' |perl -pe 's/-L $//'`
    gatk --java-options "-Xmx32g" GenomicsDBImport $gvcflist --genomicsdb-workspace-path gendb -L $interval
    gatk --java-options "-Xmx32g" GenotypeGVCFs -V gendb://gendb -R ${reffa} -D ${gatk4_dbsnp} -O gatk.vcf
    bcftools norm -c s -f ${reffa} -w 10 -O v gatk.vcf | vcf-annotate -n --fill-type gatk.vcf | bgzip > ${pair_id}.gatk.vcf.gz
    tabix ${pair_id}.gatk.vcf.gz
elif [ $algo == 'mutect' ]
then
  gatk4_dbsnp=${index_path}/clinseq_prj/dbSnp.gatk4.vcf.gz
  module load gatk/4.1.4.0
  bamlist=''
  for i in *.bam; do
      bamlist+="-I ${i} "
  done
  gatk --java-options "-Xmx20g" Mutect2 $ponopt -R ${reffa} ${bamlist} --output ${pair_id}.mutect.vcf -RF AllowAllReadsReadFilter --independent-mates  --tmp-dir `pwd`
  #gatk --java-options "-Xmx20g" FilterMutectCalls -R ${reffa} -V ${pair_id}.mutect.vcf -O ${pair_id}.mutect.filt.vcf
  vcf-sort ${pair_id}.mutect.vcf | vcf-annotate -n --fill-type | java -jar $SNPEFF_HOME/SnpSift.jar filter -p '(GEN[*].DP >= 10)' | bgzip > ${pair_id}.mutect.vcf.gz

elif [[ $algo == 'strelka2' ]]
then
    opt=''
    if [[ -n $tbed ]]
    then
	opt="--callRegions ${tbed}.gz"
    fi
    if [[ $rna == 1 ]]
    then
	mode="--rna"
    else
	mode="--exome"
    fi
    module load strelka/2.9.10 manta/1.3.1
    mkdir manta strelka
    gvcflist=''
    for i in *.bam; do
	gvcflist="$gvcflist --bam ${i}"
    done
    configManta.py $gvcflist $opt --referenceFasta ${reffa} $mode --runDir manta
    manta/runWorkflow.py -m local -j $NPROC
    if [[ -f manta/results/variants/candidateSmallIndels.vcf.gz ]]
    then
	configureStrelkaGermlineWorkflow.py $gvcflist --referenceFasta ${reffa} $mode --indelCandidates manta/results/variants/candidateSmallIndels.vcf.gz --runDir strelka
    else
	configureStrelkaGermlineWorkflow.py $gvcflist --referenceFasta ${reffa} $mode --runDir strelka
    fi
    strelka/runWorkflow.py -m local -j $NPROC
    bcftools norm -c s -f ${reffa} -w 10 -O z -o ${pair_id}.strelka2.vcf.gz strelka/results/variants/variants.vcf.gz
fi
