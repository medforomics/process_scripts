#!/bin/bash
#union.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-p  --Prefix for output file name"
  echo "Example: bash union.sh -p prefix -r /path/GRCh38"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:p:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
        h) usage;;
    esac
done
function join_by { local IFS="$1"; shift; echo "$*"; }
shift $(($OPTIND -1))
baseDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
module load bedtools/2.26.0 samtools/1.6

HS=*.hotspot.vcf.gz
list1=`ls *vcf.gz|grep -v hotspot`
list2=`ls *vcf.gz|grep -v hotspot`
varlist=''
calllist=''
for i in *.vcf.gz; do
    if [[ $i == $HS ]]
    then
	bedtools multiinter -i $list1 |cut -f 1,2,3 |bedtools intersect -header -v -a $i -b stdin |bgzip > hotspot.nooverlap.vcf.gz
	list2="$list2 hotspot.nooverlap.vcf.gz"
    fi
done 

perl $baseDir/unionvcf.pl ${index_path}/union.header.vcf $list2
perl $baseDir/vcfsorter.pl ${index_path}/genome.dict int.vcf |bgzip > ${pair_id}.union.vcf.gz

