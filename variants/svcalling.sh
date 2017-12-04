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
while getopts :r:a:b:n:k:p:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
        b) sbam=$OPTARG;;
	k) keepid=$OPTARG;;
        n) normal=$OPTARG;;
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
if [[ -z $SLURM_CPUS_ON_NODE ]]
then
    SLURM_CPUS_ON_NODE=1
fi

if [[ -a "${index_path}/genome.fa" ]]
then
    reffa="${index_path}/genome.fa"
    dict="${index_path}/genome.dict"
else 
    echo "Missing Fasta File: ${index_path}/genome.fa"
    usage

fi

module load  speedseq/20160506 novoBreak/v1.1.3 delly2/v0.7.7-multi samtools/1.6 bedtools/2.26.0 snpeff/4.3q vcftools/0.1.14
mkdir temp

if [[ -n ${normal} ]]
then
  run_novoBreak.sh /cm/shared/apps/novoBreak/novoBreak_distribution_v1.1.3rc ${reffa} ${sbam} ${normal} $SLURM_CPUS_ON_NODE
  perl $baseDir/vcf2bed.sv.pl novoBreak.pass.flt.vcf |sort -T temp -V -k 1,1 -k 2,2n > novobreak.bed
  mv novoBreak.pass.flt.vcf ${pair_id}.novobreak.vcf
  bgzip ${pair_id}.novobreak.vcf
fi
if [[ -n ${normal} ]]
then
  #RUN DELLY
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
perl $baseDir/vcf2bed.sv.pl delly.vcf > delly.bed
bgzip delly.vcf
if [[ -n ${normal} ]]
then
    tabix delly.vcf.gz
    bcftools view -O z -o delly.vcf.gz -s ${keepid} ${pair_id}.delly.vcf.gz
else
    mv delly.vcf.gz ${pair_id}.delly.vcf.gz
fi
tabix ${pair_id}.delly.vcf.gz

#MAKE FILES FOR LUMPY
samtools sort -@ $SLURM_CPUS_ON_NODE -n -o namesort.bam ${sbam}
samtools view -h namesort.bam | samblaster -M -a --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 -d discordants.sam -s splitters.sam > temp.sam
gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" splitters.sam | samtools  view -S -b - | samtools sort -o splitters.bam -
gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" discordants.sam | samtools  view -S  -b - | samtools sort -o discordants.bam -
#RUN LUMPY
if [[ -n ${normal} ]]
then
    samtools sort -@ $SLURM_CPUS_ON_NODE -n -o namesort.bam ${normal}
    samtools view -h namesort.bam | samblaster -M -a --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 -d discordants.sam -s splitters.sam > temp.sam
    gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" splitters.sam | samtools  view -S -b - | samtools sort -o normal.splitters.bam -
    gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" discordants.sam | samtools  view -S  -b - | samtools sort -o normal.discordants.bam -
    speedseq sv -t $SLURM_CPUS_ON_NODE -o sssv -R ${reffa} -B ${normal},${sbam} -D normal.discordants.bam,discordants.bam -S normal.splitters.bam,splitters.bam -x ${index_path}/exclude_alt.bed
    bcftools view -O z -o sssv.vcf.gz -s ${keepid} ${pair_id}.sssv.sv.vcf.gz 
else
    speedseq sv -t $SLURM_CPUS_ON_NODE -o ${pair_id}.sssv -R ${reffa} -B ${sbam} -D discordants.bam -S splitters.bam -x ${index_path}/exclude_alt.bed   
fi

java -jar $SNPEFF_HOME/SnpSift.jar filter "GEN[0].SU > 2" ${pair_id}.sssv.sv.vcf.gz > lumpy.vcf
perl $baseDir/vcf2bed.sv.pl lumpy.vcf > lumpy.bed

#COMPARE DELLY & LUMPY
if [[ -n ${normal} ]]
then
  bedtools multiinter -cluster -header -names novobreak delly lumpy -i novobreak.bed delly.bed lumpy.bed > sv.intersect.bed 
  grep novobreak sv.intersect.bed |cut -f 1,2,3 |sort -V -k 1,1 -k 2,2n |grep -v start | bedtools intersect -header -b stdin -a ${tid}_${nid}.novobreak.vcf.gz  | perl -p -e 's/SPIKEIN/${tid}/' |bgzip > svt1.vcf.gz
else
  bedtools multiinter -cluster -header -names delly lumpy -i delly.bed lumpy.bed > sv.intersect.bed 
fi
grep delly sv.intersect.bed |cut -f 1,2,3 |sort -V -k 1,1 -k 2,2n |grep -v 'start' |grep -v 'novobreak' | bedtools intersect -header -b stdin -a ${pair_id}.delly.vcf.gz |bgzip > svt2.vcf.gz
grep lumpy sv.intersect.bed |cut -f 1,2,3 |sort -V -k 1,1 -k 2,2n |grep -v 'start' |grep -v 'delly' |grep -v 'novobreak' | bedtools intersect -header -b stdin -a ${pair_id}.sssv.sv.vcf.gz |bgzip > svt3.vcf.gz

vcf-concat svt*.vcf.gz | vcf-sort -t temp > ${pair_id}.sv.vcf
java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config GRCh38.86 ${pair_id}.sv.vcf > ${pair_id}.sv.annot.vcf

perl $baseDir/svannot.pl -i ${pair_id}.sv.annot.vcf
bgzip ${pair_id}.sv.annot.vcf
tabix ${pair_id}.sv.annot.vcf.gz

fi
