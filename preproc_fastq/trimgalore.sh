#!/bin/bash
#trimgalore.sh

usage() {
  echo "-h Help documentation for hisat.sh"
  echo "-a  --FastQ R1"
  echo "-b  --FastQ R2"
  echo "-p  --Prefix for output file name"
  echo "Example: bash trimgalore.sh -p prefix -a SRR1551047_1.fastq.gz  -b SRR1551047_2.fastq.gz"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :a:b:p:fh opt
do
    case $opt in
        a) fq1=$OPTARG;;
        b) fq2=$OPTARG;;
        p) pair_id=$OPTARG;;
	f) filter=1;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $fq1 ]]; then
    usage
fi
fqs=("$@")
numfq=$#

if [[ -f $fq1 ]]
then
    fqs="$fq1"
    r1base="${fq1%.fastq*}"
fi
if [[ -f $fq2 ]]
then
    fqs+=" $fq2"
    r2base="${fq2%.fastq*}"
fi

numfq=${#fqs[@]}

source /etc/profile.d/modules.sh
module load trimgalore/0.6.4 cutadapt/2.5 perl/5.28.0 

if [ $numfq > 1 ]
then
    trim_galore --paired -q 25 --illumina --gzip --length 35 ${fqs}
    mv ${r1base}_val_1.fq.gz ${pair_id}.trim.R1.fastq.gz
    mv ${r2base}_val_2.fq.gz ${pair_id}.trim.R2.fastq.gz
else
    trim_galore -q 25 --illumina --gzip --length 35 ${fqs}
    mv ${r1base}_trimmed.fq.gz ${pair_id}.trim.R1.fastq.gz
    cp ${pair_id}.trim.R1.fastq.gz ${pair_id}.trim.R2.fastq.gz 
fi

if [[ $filter == 1 ]]
then
      perl $baseDir/parse_trimreport.pl ${pair_id}.trimreport.txt *trimming_report.txt
fi
