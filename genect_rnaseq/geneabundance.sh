#!/bin/bash
#rnaseqalign.sh

usage() {
  echo "-h Help documentation for hisat.sh"
  echo "-b  --BAM File"
  echo "-g  --GTF File"
  echo "-s  --stranded"
  echo "-p  --Prefix for output file name"
  echo "Example: bash geneabudance.sh genect_rnaseq/geneabundance.sh -s stranded -g gtf_file -p pair_id -b bam"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :b:g:p:f:s:h opt
do
    case $opt in
        b) sbam=$OPTARG;;
        g) gtf=$OPTARG;;
        p) pair_id=$OPTARG;;
        s) stranded=$OPTARG;;
	f) filter=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $sbam ]]
then
    usage
fi
NPROC=$SLURM_CPUS_ON_NODE
if [[ -z $NPROC ]]
then
    NPROC=`nproc`
fi
if [[ $NPROC > 64 ]]
then
    NPROC=64
fi
source /etc/profile.d/modules.sh
module load subread/1.6.1

export PATH=/project/shared/bicf_workflow_ref/seqprg/bin:$PATH

featureCounts -s $stranded -M --fraction -J --ignoreDup -T $NPROC -p -g gene_name -a ${gtf} -o ${pair_id}.cts ${sbam}

mkdir ${pair_id}_stringtie
cd ${pair_id}_stringtie
stringtie ../${sbam} -p $NPROC -G ${gtf} -B -e -o denovo.gtf -A ../${pair_id}.fpkm.txt

if [[ -f $filter ]]
then
    cd ..
    mv ${pair_id}.fpkm.txt ${prefix}.fpkm.ori.txt
    perl ${baseDir}/fpkm_subset_panel.pl -f ${prefix}.fpkm.ori.txt -g $filter
    mv ${prefix}.fpkm.capture.txt ${pair_id}.fpkm.txt
fi
