#!/bin/bash
#annotvcf.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-p  --Prefix for output file name"
  echo "-v  --VCF File"
  echo "Example: bash hisat.sh -p prefix -r /path/GRCh38 -v vcffile"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:v:p:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
	v) unionvcf=$OPTARG;;
        h) usage;;
    esac
done
function join_by { local IFS="$1"; shift; echo "$*"; }
shift $(($OPTIND -1))

source /etc/profile.d/modules.sh
module load python/2.7.x-anaconda bedtools/2.26.0 samtools/1.6 snpeff/4.3q

if [[ $index_path == '/project/shared/bicf_workflow_ref/GRCh38/hisat_index' ]]
then
    index_path='/project/shared/bicf_workflow_ref/GRCh38'
fi

if  [[ $index_path == '/project/shared/bicf_workflow_ref/GRCh38' ]] 
then
    tabix ${unionvcf}
    bcftools annotate -Oz -a ${index_path}/ExAC.vcf.gz -o ${pair_id}.exac.vcf.gz --columns CHROM,POS,AC_Het,AC_Hom,AC_Hemi,AC_Adj,AN_Adj,AC_POPMAX,AN_POPMAX,POPMAX ${unionvcf}
    tabix ${pair_id}.exac.vcf.gz 
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config GRCh38.86 ${pair_id}.exac.vcf.gz | java -jar $SNPEFF_HOME/SnpSift.jar annotate -id ${index_path}/dbSnp.vcf.gz -  | java -jar $SNPEFF_HOME/SnpSift.jar annotate -info CLNSIG,CLNDSDB,CLNDSDBID,CLNDBN,CLNREVSTAT,CLNACC ${index_path}/clinvar.vcf.gz - | java -jar $SNPEFF_HOME/SnpSift.jar annotate -info CNT ${index_path}/cosmic.vcf.gz - | java -Xmx10g -jar $SNPEFF_HOME/SnpSift.jar dbnsfp -v -db ${index_path}/dbNSFP.txt.gz - |bgzip > ${pair_id}.annot.vcf.gz
    tabix ${pair_id}.annot.vcf.gz
else 
    if [[ $index_path == '/project/shared/bicf_workflow_ref/GRCm38' ]]
    then
	snpeffvers='GRCh38.86'
    elif [[ $index_path == '/project/shared/bicf_workflow_ref/GRCh37' ]]
    then
	snpeffvers='GRCh37.75'
    fi
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffvers} ${unionvcf} |bgzip > ${pair_id}.annot.vcf.gz
    tabix ${pair_id}.annot.vcf.gz
fi
