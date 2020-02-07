#!/bin/bash
#svcalling.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Path to Reference Genome with the file genome.fa"
  echo "-p  --Prefix for output file name"
  echo "-b  --Bam File"
  echo "-n  --Reference Bam File"
  echo "Example: bash svcalling.sh -p prefix -r /path/GRCh38 -a gatk"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:p:b:i:x:y:n:l:a:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
        b) sbam=$OPTARG;;
        i) tumor=$OPTARG;;
        n) normal=$OPTARG;;
	a) method=$OPTARG;;
        x) tid=$OPTARG;;
        y) nid=$OPTARG;;
	l) bed=$OPTARG;;
        h) usage;;
    esac
done
function join_by { local IFS="$1"; shift; echo "$*"; }

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"


# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $index_path ]]; then
    usage
fi
NPROC=$SLURM_CPUS_ON_NODE
if [[ -z $NPROC ]]
then
    NPROC=`nproc`
fi

if [[ -a "${index_path}/genome.fa" ]]
then
    reffa="${index_path}/genome.fa"
    dict="${index_path}/genome.dict"
else 
    echo "Missing Fasta File: ${index_path}/genome.fa"
    usage

fi
if [[ -z $snpeffgeno ]]
then
    snpeffgeno='GRCh38.86'
fi

source /etc/profile.d/modules.sh	
module load htslib/gcc/1.8 samtools/gcc/1.8 bcftools/gcc/1.8 bedtools/2.26.0 snpeff/4.3q vcftools/0.1.14

if [[ $method == 'delly' ]]
then
    mkdir temp
    module load  delly2/v0.7.7-multi
    if [[ -n ${normal} ]]
    then
	#RUN DELLY
	echo -e "${nid}\tcontrol"> samples.tsv
	echo -e "${tid}\ttumor" >> samples.tsv
	delly2 call -t BND -o delly_translocations.bcf -q 30 -g ${reffa} ${sbam} ${normal}
	delly2 call -t DUP -o delly_duplications.bcf -q 30 -g ${reffa} ${sbam} ${normal}
	delly2 call -t INV -o delly_inversions.bcf -q 30 -g ${reffa} ${sbam} ${normal}
	delly2 call -t DEL -o delly_deletion.bcf -q 30 -g ${reffa} ${sbam} ${normal}
	delly2 call -t INS -o delly_insertion.bcf -q 30 -g ${reffa} ${sbam} ${normal}
	delly2 filter -t BND -o  delly_tra.bcf -f somatic -s samples.tsv delly_translocations.bcf
	delly2 filter -t DUP -o  delly_dup.bcf -f somatic -s samples.tsv delly_duplications.bcf
	delly2 filter -t INV -o  delly_inv.bcf -f somatic -s samples.tsv delly_inversions.bcf
	delly2 filter -t DEL -o  delly_del.bcf -f somatic -s samples.tsv delly_deletion.bcf
	delly2 filter -t INS -o  delly_ins.bcf -f somatic -s samples.tsv delly_insertion.bcf
    else
	#RUN DELLY
	delly2 call -t BND -o delly_translocations.bcf -q 30 -g ${reffa} ${sbam}
	delly2 call -t DUP -o delly_duplications.bcf -q 30 -g ${reffa} ${sbam}
	delly2 call -t INV -o delly_inversions.bcf -q 30 -g ${reffa} ${sbam}
	delly2 call -t DEL -o delly_deletion.bcf -q 30 -g ${reffa} ${sbam}
	delly2 call -t INS -o delly_insertion.bcf -q 30 -g ${reffa} ${sbam}
	delly2 filter -t BND -o  delly_tra.bcf -f germline delly_translocations.bcf
	delly2 filter -t DUP -o  delly_dup.bcf -f germline delly_duplications.bcf
	delly2 filter -t INV -o  delly_inv.bcf -f germline delly_inversions.bcf
	delly2 filter -t DEL -o  delly_del.bcf -f germline delly_deletion.bcf
	delly2 filter -t INS -o  delly_ins.bcf -f germline delly_insertion.bcf
    fi
    #MERGE DELLY AND MAKE BED
    bcftools concat -a -O v delly_dup.bcf delly_inv.bcf delly_tra.bcf delly_del.bcf delly_ins.bcf| vcf-sort -t temp > delly.vcf
    bgzip delly.vcf
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} delly.vcf | bgzip > ${pair_id}.delly.vcf.gz
fi
if [[ $method == 'svaba' ]]
then
    if [[ -n ${normal} ]]
    then
	/project/shared/bicf_workflow_ref/seqprg/svaba/bin/svaba run -p $NPROC -G ${reffa} -t ${sbam} -n ${normal} -a ${pair_id}
    else
	/project/shared/bicf_workflow_ref/seqprg/svaba/bin/svaba run -p $NPROC -G ${reffa} -t ${sbam} -a ${pair_id}
    fi
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} ${pair_id}.svaba.unfiltered.somatic.sv.vcf | bgzip > ${pair_id}.svaba.vcf.gz
fi

if [[ $method == 'lumpy' ]]
then
    #MAKE FILES FOR LUMPY
    samtools sort -@ $NPROC -n -o namesort.bam ${sbam}
    samtools view -h namesort.bam | samblaster -M -a --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 -d discordants.sam -s splitters.sam > temp.sam
    gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" splitters.sam | samtools  view -S -b - | samtools sort -o splitters.bam -
    gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" discordants.sam | samtools  view -S  -b - | samtools sort -o discordants.bam -
    #RUN LUMPY
    if [[ -n ${normal} ]]
    then
	samtools sort -@ $NPROC -n -o namesort.bam ${normal}
	samtools view -h namesort.bam | samblaster -M -a --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 -d discordants.sam -s splitters.sam > temp.sam
	gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" splitters.sam | samtools  view -S -b - | samtools sort -o normal.splitters.bam -
	gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" discordants.sam | samtools  view -S  -b - | samtools sort -o normal.discordants.bam -
	speedseq sv -t $NPROC -o lumpy -R ${reffa} -B ${normal},${sbam} -D normal.discordants.bam,discordants.bam -S normal.splitters.bam,splitters.bam -x ${index_path}/exclude_alt.bed
    else
	speedseq sv -t $NPROC -o lumpy -R ${reffa} -B ${sbam} -D discordants.bam -S splitters.bam -x ${index_path}/exclude_alt.bed   
    fi
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} lumpy.sv.vcf.gz | java -jar $SNPEFF_HOME/SnpSift.jar filter " ( GEN[*].DV >= 20 )" | bgzip > ${pair_id}.lumpy.vcf.gz
fi
if [[ $method == 'pindel' ]]
then
    module load pindel/0.2.5-intel
    genomefiledate=`find ${reffa} -maxdepth 0 -printf "%TY%Tm%Td\n"`
    touch ${pair_id}.pindel.config
    for i in *.bam; do
	sname=`echo "$i" |cut -f 1 -d '.'`
	echo -e "${i}\t400\t${sname}" >> ${pair_id}.pindel.config
	samtools index -@ $NPROC $i
    done
    pindel -T $NPROC -f ${reffa} -i ${pair_id}.pindel.config -o ${pair_id}.pindel_out --RP
    pindel2vcf -P ${pair_id}.pindel_out -r ${reffa} -R HG38 -d ${genomefiledate} -v pindel.vcf
    cat pindel.vcf | java -jar $SNPEFF_HOME/SnpSift.jar filter "( GEN[*].AD[1] >= 10 )" | bgzip > pindel.vcf.gz
    tabix pindel.vcf.gz
    bash $baseDir/norm_annot.sh -r ${index_path} -p pindel -v pindel.vcf.gz
    perl $baseDir/parse_pindel.pl ${pair_id} pindel.norm.vcf.gz
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} ${pair_id}.indel.vcf |bgzip > ${pair_id}.pindel_indel.vcf.gz
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} ${pair_id}.dup.vcf | bedtools intersect -header -b ${bed} -a stdin | bgzip > ${pair_id}.pindel_tandemdup.vcf.gz
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} ${pair_id}.sv.vcf | bgzip > ${pair_id}.pindel_sv.vcf.gz
fi
if [[ $method == 'itdseek' ]]
then
    stexe=`which samtools`
    samtools view -@ $NPROC -L ${bed} ${sbam} | /project/shared/bicf_workflow_ref/seqprg/itdseek-1.2/itdseek.pl --refseq ${reffa} --samtools ${stexe} --bam ${sbam} | vcf-sort | bedtools intersect -header -b ${bed} -a stdin | bgzip > ${pair_id}.itdseek.vcf.gz

    tabix ${pair_id}.itdseek.vcf.gz
    bcftools norm --fasta-ref $reffa -m - -Ov ${pair_id}.itdseek.vcf.gz | java -Xmx30g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} - |bgzip > ${pair_id}.itdseek_tandemdup.vcf.gz
fi
