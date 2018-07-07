#!/bin/bash
#union.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-p  --Prefix for output file name"
  echo "-v  --VCF File" 
  echo "Example: bash union.sh -p prefix -r /path/GRCh38"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:p:v:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
	v) vcf=$OPTARG;;
        h) usage;;
    esac
done
function join_by { local IFS="$1"; shift; echo "$*"; }
shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

source /etc/profile.d/modules.sh
module load bedtools/2.26.0 samtools/1.6 bcftools/1.6 snpeff/4.3q 

prefix="${vcf%.vcf.gz}"
perl $baseDir\/uniform_vcf_gt.pl $vcf
bgzip -f ${prefix}.uniform.vcf
j=${prefix}.uniform.vcf.gz
tabix -f $j
bcftools norm -m - -O z -o ${prefix}.norm.vcf.gz $j
