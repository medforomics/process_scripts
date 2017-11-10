!/bin/bash
#svcalling.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Path to Reference Genome with the file genome.fa"
  echo "-p  --Prefix for output file name"
  echo "-b  --Bam File"
  echo "Example: bash svcalling.sh -p prefix -r /path/GRCh38 -a gatk"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:a:b:p:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
        b) sbam=$OPTARG;;
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
baseDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

module load  speedseq/20160506 novoBreak/v1.1.3 delly2/v0.7.7-multi samtools/1.6 bedtools/2.26.0 snpeff/4.3q vcftools/0.1.14
#RUN DELLY
mkdir temp
delly2 call -t BND -o delly_translocations.bcf -q 30 -g ${reffa} ${sbam}
delly2 call -t DUP -o delly_duplications.bcf -q 30 -g ${reffa} ${sbam}
delly2 call -t INV -o delly_inversions.bcf -q 30 -g ${reffa} ${sbam}
delly2 call -t DEL -o delly_deletion.bcf -q 30 -g ${reffa} ${sbam}
delly2 call -t INS -o delly_insertion.bcf -q 30 -g ${reffa} ${sbam}
delly2 filter -t BND -o  delly_tra.bcf -f germline delly_translocations.bcf
delly2 filter -t DUP -o  delly_dup.bcf -f germline delly_translocations.bcf
delly2 filter -t INV -o  delly_inv.bcf -f germline delly_translocations.bcf
delly2 filter -t DEL -o  delly_del.bcf -f germline delly_translocations.bcf
delly2 filter -t INS -o  delly_ins.bcf -f germline delly_translocations.bcf

#MERGE DELLY AND MAKE BED
bcftools concat -a -O v delly_dup.bcf delly_inv.bcf delly_tra.bcf delly_del.bcf delly_ins.bcf | vcf-sort > ${pair_id}.delly.vcf
perl $baseDir/vcf2bed.sv.pl ${pair_id}.delly.vcf > delly.bed
bgzip ${pair_id}.delly.vcf
tabix ${pair_id}.delly.vcf.gz

#MAKE FILES FOR LUMPY
samtools sort -@ $SLURM_CPUS_ON_NODE -n -o namesort.bam ${sbam}
samtools view -h namesort.bam | samblaster -M -a --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 -d discordants.sam -s splitters.sam > temp.sam
gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" splitters.sam | samtools  view -S -b - | samtools sort -o splitters.bam -
gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" discordants.sam | samtools  view -S  -b - | samtools sort -o discordants.bam -

#RUN LUMPY
speedseq sv -t $SLURM_CPUS_ON_NODE -o ${pair_id}.sssv -R ${reffa} -B ${sbam} -D discordants.bam -S splitters.bam -x ${index_path}/exclude_alt.bed
java -jar $SNPEFF_HOME/SnpSift.jar filter "GEN[0].SU > 2" ${pair_id}.sssv.sv.vcf.gz > lumpy.vcf
perl $baseDir/vcf2bed.sv.pl lumpy.vcf > lumpy.bed

#COMPARE DELLY & LUMPY
bedtools intersect -v -a lumpy.bed -b delly.bed > lumpy_only.bed
bedtools intersect -header -b lumpy_only.bed -a lumpy.vcf |bgzip > lumpy_only.vcf.gz
vcf-concat ${pair_id}.delly.vcf.gz lumpy_only.vcf.gz |vcf-sort -t temp > ${pair_id}.sv.vcf
perl $baseDir/vcf2bed.sv.pl ${pair_id}.sv.vcf |sort -V -k 1,1 -k 2,2n | grep -v 'alt' |grep -v 'random' |uniq > svs.bed
bedtools intersect -header -wb -a svs.bed -b ${index_path}/gencode.exons.bed > exonoverlap_sv.txt
bedtools intersect -v -header -wb -a svs.bed -b ${index_path}/gencode.exons.bed | bedtools intersect -header -wb -a stdin -b ${index_path}/gencode.genes.chr.bed > geneoverlap_sv.txt
perl $baseDir/annot_sv.pl -r ${index_path} -i ${pair_id}.sv.vcf
bgzip ${pair_id}.sv.vcf

