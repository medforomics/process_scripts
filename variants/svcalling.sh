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
while getopts :r:p:b:t:x:c:g:y:n:l:a:hf opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
        t) tumor=$OPTARG;;
        n) normal=$OPTARG;;
	a) method=$OPTARG;;
        x) tid=$OPTARG;;
        y) nid=$OPTARG;;
	f) filter=1;;
	g) snpeffgeno=$OPTARG;;
        b) sbam=$OPTARG;;
	c) tbed=$OPTARG;;
	l) itdbed=$OPTARG;;
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
if [[ -z $isdocker ]]
then
    source /etc/profile.d/modules.sh	
    module load htslib/gcc/1.8 samtools/gcc/1.8 bcftools/gcc/1.8 bedtools/2.26.0 snpeff/4.3q vcftools/0.1.14
    export PATH=/project/shared/bicf_workflow_ref/seqprg/bin:$PATH
fi
mkdir -p temp

if [[ -z $tid ]]
then
    tid=`samtools view -H ${sbam} |grep '^@RG' |perl -pe 's/\t/\n/g' |grep ID |cut -f 2 -d ':'`
fi

bams=''
for i in *.bam; do
    bams="$bams $i" 
done

#RUN DELLY
if [[ $method == 'delly' ]]
then
    delly2 call -t BND -o ${pair_id}.delly_translocations.bcf -q 30 -g ${reffa} ${bams} 
    delly2 call -t DUP -o ${pair_id}.delly_duplications.bcf -q 30 -g ${reffa} ${bams}
    delly2 call -t INV -o ${pair_id}.delly_inversions.bcf -q 30 -g ${reffa} ${bams}
    delly2 call -t DEL -o ${pair_id}.delly_deletion.bcf -q 30 -g ${reffa} ${bams}
    delly2 call -t INS -o ${pair_id}.delly_insertion.bcf -q 30 -g ${reffa} ${bams}
    #MERGE DELLY AND MAKE BED
    bcftools concat -a -O v ${pair_id}.delly_duplications.bcf ${pair_id}.delly_inversions.bcf ${pair_id}.delly_translocations.bcf ${pair_id}.delly_deletion.bcf ${pair_id}.delly_insertion.bcf | vcf-sort -t temp | bgzip > ${pair_id}.delly.svar.vcf.gz
    bash $baseDir/norm_annot.sh -r ${index_path} -p ${pair_id}.delly.sv -v ${pair_id}.delly.svar.vcf.gz -s
    java -jar $SNPEFF_HOME/SnpSift.jar filter "( GEN[*].DP >= 20 )" ${pair_id}.delly.sv.norm.vcf.gz | java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} - | bgzip > ${pair_id}.delly.vcf.gz
    if [[ $filter == 1 ]]
    then
	zcat ${pair_id}.delly.vcf.gz | $SNPEFF_HOME/scripts/vcfEffOnePerLine.pl |java -jar $SNPEFF_HOME/SnpSift.jar extractFields - CHROM POS CHR2 END ANN[*].EFFECT ANN[*].GENE ANN[*].BIOTYPE FILTER FORMAT GEN[*] |grep -E 'gene_fusion|feature_fusion' | sort -u > ${pair_id}.dgf.txt
	mv ${pair_id}.delly.vcf.gz ${pair_id}.delly.ori.vcf.gz
	perl $baseDir/filter_delly.pl -t $tid -p $pair_id -i ${pair_id}.delly.ori.vcf.gz
	bgzip -f ${pair_id}.delly.vcf
	zgrep '#CHROM' ${pair_id}.delly.vcf.gz > ${pair_id}.delly.genefusion.txt
	cat ${pair_id}.delly.potentialfusion.txt  ${pair_id}.dgf.txt |sort -u >> ${pair_id}.delly.genefusion.txt
    fi
elif [[ $method == 'svaba' ]]
then
    if [[ -n ${normal} ]]
    then
	svaba run -p $NPROC -G ${reffa} -t ${sbam} -n ${normal} -a ${pair_id}
    else
	bamlist=''
	for i in *.bam; do
	    bamlist="$bamlist -t ${i}"
	done
	svaba run -p $NPROC -G ${reffa} -a ${pair_id} $bamlist
    fi
    #Create SV FILE
    vcf-concat ${pair_id}.svaba.unfiltered*sv.vcf | perl -pe 's/\.consensus|\.bam//g' | vcf-sort| bgzip > ${pair_id}.svaba.unfiltered.sv.vcf.gz
    bash $baseDir/norm_annot.sh -r ${index_path} -p svaba.sv -v ${pair_id}.svaba.unfiltered.sv.vcf.gz
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} svaba.sv.norm.vcf.gz | java -jar $SNPEFF_HOME/SnpSift.jar filter "( GEN[*].AO >= 20)" | bgzip > ${pair_id}.svaba.sv.vcf.gz

    vcf-concat ${pair_id}.svaba.unfiltered*indel.vcf | perl -pe 's/\.consensus|\.bam//g' | vcf-sort | java -jar $SNPEFF_HOME/SnpSift.jar filter "( SPAN >= 20)" - |bgzip > ${pair_id}.svaba.indel.vcf.gz
    bash $baseDir/norm_annot.sh -r ${index_path} -p svaba.indel -v ${pair_id}.svaba.indel.vcf.gz
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} svaba.indel.norm.vcf.gz | bgzip > ${pair_id}.svaba.vcf.gz

    if [[ $filter == 1 ]]
    then
	zcat ${pair_id}.svaba.sv.vcf.gz | $SNPEFF_HOME/scripts/vcfEffOnePerLine.pl |java -jar $SNPEFF_HOME/SnpSift.jar extractFields - CHROM POS ALT ID ANN[*].EFFECT ANN[*].GENE ANN[*].BIOTYPE FILTER FORMAT GEN[*] |grep -E 'gene_fusion|feature_fusion' | sort -u  > ${pair_id}.sgf.txt
	mv ${pair_id}.svaba.vcf.gz ${pair_id}.svaba.ori.vcf.gz
	perl $baseDir/filter_svaba.pl -t $tid -p ${pair_id} -i ${pair_id}.svaba.ori.vcf.gz -s ${pair_id}.svaba.sv.vcf.gz
	bgzip ${pair_id}.svaba.vcf
	zgrep '#CHROM' ${pair_id}.svaba.sv.vcf.gz > ${pair_id}.svaba.genefusion.txt
	cat ${pair_id}.svaba.potentialfusion.txt ${pair_id}.sgf.txt | sort -u >> ${pair_id}.svaba.genefusion.txt
       fi
elif [[ $method == 'lumpy' ]]
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
elif [[ $method == 'pindel' ]]
then
    if [[ -z $isdocker ]]
    then
	module load pindel/0.2.5-intel
    fi
    genomefiledate=`find ${reffa} -maxdepth 0 -printf "%TY%Tm%Td\n"`
    touch ${pair_id}.pindel.config
    for i in *.bam; do
	sname=`echo "$i" |cut -f 1 -d '.'`
	echo -e "${i}\t400\t${sname}" >> ${pair_id}.pindel.config
	samtools index -@ $NPROC $i
    done
    bedopt=''
    if [[ -f $tbed ]]
    then
	bedopt="-j $tbed"
    fi
    pindel -T $NPROC -f ${reffa} -i ${pair_id}.pindel.config -o ${pair_id}.pindel_out --report_inversions false $bedopt
    pindel2vcf -P ${pair_id}.pindel_out -r ${reffa} -R HG38 -d ${genomefiledate} -v pindel.vcf
    cat pindel.vcf | java -jar $SNPEFF_HOME/SnpSift.jar filter "( GEN[*].AD[1] >= 10 )" | bgzip > pindel.vcf.gz
    tabix pindel.vcf.gz
    bash $baseDir/norm_annot.sh -r ${index_path} -p pindel -v pindel.vcf.gz
    perl $baseDir/parse_pindel.pl ${pair_id} pindel.norm.vcf.gz
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} ${pair_id}.indel.vcf |bgzip > ${pair_id}.pindel_indel.vcf.gz
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} ${pair_id}.dup.vcf | bedtools intersect -header -b ${itdbed} -a stdin | bgzip > ${pair_id}.pindel_tandemdup.vcf.gz
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} ${pair_id}.sv.vcf | bgzip > ${pair_id}.pindel.sv.vcf.gz
    if [[ $filter == 1 ]]
    then
	perl $baseDir/filter_pindel.pl -d ${pair_id}.pindel_tandemdup.vcf.gz -s ${pair_id}.pindel.sv.vcf.gz -i ${pair_id}.pindel_indel.vcf.gz
	mv ${pair_id}.pindel_tandemdup.vcf.gz ${pair_id}.pindel_tandemdup.unfilt.vcf.gz
	mv ${pair_id}.pindel_tandemdup.pass.vcf ${pair_id}.pindel_tandemdup.vcf
	bgzip ${pair_id}.pindel_tandemdup.vcf
	mv ${pair_id}.pindel_indel.pass.vcf ${pair_id}.pindel.vcf
	bgzip ${pair_id}.pindel.vcf
	mv ${pair_id}.pindel.sv.vcf.gz ${pair_id}.pindel.sv.unfilt.vcf.gz
	mv ${pair_id}.pindel.sv.pass.vcf ${pair_id}.pindel.sv.vcf
	bgzip ${pair_id}.pindel.sv.vcf
	zgrep '#CHROM' ${pair_id}.pindel.sv.vcf.gz > ${pair_id}.pindel.genefusion.txt
	zcat ${pair_id}.pindel.sv.vcf.gz | $SNPEFF_HOME/scripts/vcfEffOnePerLine.pl |java -jar $SNPEFF_HOME/SnpSift.jar extractFields - CHROM POS CHROM END ANN[*].EFFECT ANN[*].GENE ANN[*].BIOTYPE  FILTER FORMAT GEN[*] |grep -E 'gene_fusion|feature_fusion' | sort -u >> ${pair_id}.pindel.genefusion.txt
    fi
elif [[ $method == 'itdseek' ]]
then
    stexe=`which samtools`
    echo $stexe
    samtools view -@ $NPROC -L ${itdbed} ${sbam} | itdseek.pl --refseq ${reffa} --samtools ${stexe} --bam ${sbam} | vcf-sort | bedtools intersect -header -b ${itdbed} -a stdin | java -Xmx30g -jar $SNPEFF_HOME/SnpSift.jar filter "( LEN < 10000 )" | bgzip > ${pair_id}.itdseek.vcf.gz
    
    tabix ${pair_id}.itdseek.vcf.gz
    
    bcftools norm --fasta-ref $reffa -m - -Ov ${pair_id}.itdseek.vcf.gz | java -Xmx30g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} - |bgzip > ${pair_id}.itdseek_tandemdup.vcf.gz
    if [[ $filter == 1 ]]
    then
	perl $baseDir/filter_itdseeker.pl -t ${pair_id} -d ${pair_id}.itdseek_tandemdup.vcf.gz
	mv ${pair_id}.itdseek_tandemdup.vcf.gz ${pair_id}.itdseek_tandemdup.unfilt.vcf.gz
	mv ${pair_id}.itdseek_tandemdup.pass.vcf ${pair_id}.itdseek_tandemdup.vcf
	bgzip ${pair_id}.itdseek_tandemdup.vcf
    fi
fi
