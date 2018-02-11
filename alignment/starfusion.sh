#!/bin/bash
#trimgalore.sh

usage() {
  echo "-h Help documentation for hisat.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-a  --FastQ R1"
  echo "-b  --FastQ R2"
  echo "-p  --Prefix for output file name"
  echo "Example: bash starfusion.sh -p prefix -r /project/shared/bicf_workflow_ref/GRCh38 -a SRR1551047_1.fastq.gz  -b SRR1551047_2.fastq.gz"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:a:b:p:fh opt
do
    case $opt in
        r) refgeno=$OPTARG;;
        a) fq1=$OPTARG;;
        b) fq2=$OPTARG;;
        p) pair_id=$OPTARG;;
	f) filter=1;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $fq1 ]]; then
    usage
fi

index_path=${refgeno}/CTAT_lib/
baseDir="`dirname \"$0\"`"
source /etc/profile.d/modules.sh
module add python/2.7.x-anaconda star/2.5.2b
STAR-Fusion --genome_lib_dir ${index_path} --left_fq ${fq1} --right_fq ${fq2} --output_dir star_fusion &> star_fusion.err
mv star_fusion/star-fusion.fusion_candidates.final.abridged ${pair_id}.starfusion.txt

if [[ $filter==1 ]]
then
perl $baseDir/filter_genefusions.pl -p ${pair_id} -f ${pair_id}.starfusion.txt
fi
