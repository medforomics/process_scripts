#!/bin/bash
#cnvkit.sh

usage() {
  echo "-h  --Help documentation for markdups.sh"
  echo "-b  --BAM file"
  echo "-p  --Prefix for output file name"
  echo "Example: bash cnvkit.sh -p prefix -b file.bam"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :b:p:uh opt
do
    case $opt in
        b) sbam=$OPTARG;;
        p) pair_id=$OPTARG;;
	u) umi='umi';;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $sbam ]]; then
    usage
fi
if [[ -z $SLURM_CPUS_ON_NODE ]]
then
    SLURM_CPUS_ON_NODE=1
fi
baseDir="`dirname \"$0\"`"
module load cnvkit/0.9.0
cnvkit.py coverage ${sbam} /project/shared/bicf_workflow_ref/GRCh38/UTSWV2.cnvkit_targets.bed -o ${pair_id}.targetcoverage.cnn
cnvkit.py coverage ${sbam} /project/shared/bicf_workflow_ref/GRCh38/UTSWV2.cnvkit_antitargets.bed -o ${pair_id}.antitargetcoverage.cnn
if [[ $umi == 'umi' ]]
then
cnvkit.py fix ${pair_id}.targetcoverage.cnn ${pair_id}.antitargetcoverage.cnn /project/shared/bicf_workflow_ref/GRCh38/UTSWV2.uminormals.cnn -o ${pair_id}.cnr
else
cnvkit.py fix ${pair_id}.targetcoverage.cnn ${pair_id}.antitargetcoverage.cnn /project/shared/bicf_workflow_ref/GRCh38/UTSWV2.normals.cnn -o ${pair_id}.cnr
fi   

cnvkit.py segment ${pair_id}.cnr -o ${pair_id}.cns
cnvkit.py call ${pair_id}.cns -o ${pair_id}.call.cns
cnvkit.py diagram ${pair_id}.cnr -s ${pair_id}.cns -o ${pair_id}.cnv.pdf
perl $baseDir/filter_cnvkit.pl *.call.cns
