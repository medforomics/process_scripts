#!/bin/bash
#cnvkit.sh

usage() {
  echo "-h  --Help documentation for markdups.sh"
  echo "-b  --BAM file"
  echo "-p  --Prefix for output file name"
  echo "-n  --Panel of Normal cnn file"
  echo "-t  --Target and Antitarget prefix"
  echo "Example: bash cnvkit.sh -p prefix -b file.bam"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :b:p:n:t:r:uqh opt
do
    case $opt in
        b) sbam=$OPTARG;;
        p) pair_id=$OPTARG;;
	n) normals=$OPTARG;;
	r) index_path=$OPTARG;;
	t) targets=$OPTARG;;
	u) umi='umi';;
	q) idtsnp=1;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

if [[ -z $index_path ]] 
then
    index_path='/project/shared/bicf_workflow_ref/human/grch38_cloud/dnaref'
fi
# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $sbam ]]
then
    usage
fi
if [[ -z $SLURM_CPUS_ON_NODE ]]
then
    SLURM_CPUS_ON_NODE=1
fi
if [[ -z $normals ]] || [[ -z $targets ]]
then
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

echo "${targets}targets.bed"
echo "${targets}antitargets.bed"
echo "${normals}"

source /etc/profile.d/modules.sh
module load cnvkit/0.9.5 bedtools/2.26.0 samtools/gcc/1.8 bcftools/gcc/1.8
unset DISPLAY

if [[ $idtsnp == 1 ]]
then
    samtools index ${sbam}
    bcftools mpileup --threads 10 --gvcf 10 -A -a 'INFO/AD,INFO/ADF,INFO/ADR,FORMAT/DP,FORMAT/SP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR' -Ou -Q20 -d 1000000 -L 1000000 -C50 -f ${reffa} ${sbam} | bcftools call --threads 10 -vmO v -o common_variants.vcf -T ${index_path}/IDT_snps.hg38.bed
    $baseDir/formatVcfCNV.pl cnvkit_common common_variants.vcf
fi

cnvkit.py coverage ${sbam} ${targets}targets.bed -o ${pair_id}.targetcoverage.cnn
cnvkit.py coverage ${sbam} ${targets}antitargets.bed -o ${pair_id}.antitargetcoverage.cnn
cnvkit.py fix ${pair_id}.targetcoverage.cnn ${pair_id}.antitargetcoverage.cnn ${normals} -o ${pair_id}.cnr
cnvkit.py segment ${pair_id}.cnr -o ${pair_id}.cns

if [[ $idtsnp == 1 ]]
then
    cnvkit.py call --filter cn ${pair_id}.cns -v cnvkit_common.vcf -o ${pair_id}.call.cns
else 
    cnvkit.py call --filter cn ${pair_id}.cns -o ${pair_id}.call.cns
fi

cnvkit.py scatter ${pair_id}.cnr -s ${pair_id}.call.cns -t --segment-color "blue" -o ${pair_id}.cnv.scatter.pdf
cut -f 1,2,3 ${pair_id}.call.cns | grep -v chrom | bedtools intersect -wao -b /project/shared/bicf_workflow_ref/human/GRCh38/cytoBand.txt -a stdin |cut -f 1,2,3,7 >  ${pair_id}.cytoband.bed
perl $baseDir/filter_cnvkit.pl ${pair_id}.call.cns
