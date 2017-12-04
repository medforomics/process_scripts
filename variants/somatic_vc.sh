#!/bin/bash
#run_somatic_caller.sh


usage(){
  echo "-h --Help documentation for run_somatic_caller.sh"
  echo "-a --Somatic Workflow Method: strelka2, virmid, speedseq, mutect2, varscan, shimmer, lancet"
  echo "-r  --Path to Reference Genome with the file genome.fa"
  echo "-n --Normal"
  echo "-t --Tumor"
  echo "-x --NormalID"
  echo "-y --TumorID"
  echo "-i --NormalBAM used for Mantra in the case of UMI consensus"
  echo "-j --TumorBAM used for Mantra in the case of UMI consensus"
  echo "Example: bash somatic_vc.sh -a strelka2 -y ORD1_N_panel1385 -y ORD1_T_panel138 -n ORD1_N_panel1385.final.bam -t ORD1_T_panel1385.final.bam"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :n:t:r:x:y:i:j:a:h opt
do
    case $opt in 
      r) index_path=$OPTARG;;
      x) tid=$OPTARG;;
      y) nid=$OPTARG;;
      n) normal=$OPTARG;;
      t) tumor=$OPTARG;;
      i) mnormal=$OPTARG;;
      j) mtumor=$OPTARG;;
      a) algo=$OPTARG;;
      h) usage;;
    esac
done

shift $(($OPTIND -1))

#Check for mandatory options
if [[ -z $normal ]] || [[ -z $tumor ]] || [[ -z $algo ]]; then
    echo $normal $tumor $algo
    usage
fi 
if [[ -z $SLURM_CPUS_ON_NODE ]] 
  then 
    SLURM_CPUS_ON_NODE=1
fi
pair_id=${tid}_${nid}
if [[ -z $mtumor ]]
then
    mtumor=tumor
    mnormal=normal
fi

if [[ -a "${index_path}/genome.fa" ]]
then
    reffa="${index_path}/genome.fa"
    dict="${index_path}/genome.dict"
else 
    echo "Missing Fasta File: ${index_path}/genome.fa"
    usage
fi
if [[ -a "${index_path}/dbSnp.vcf.gz" ]]
then
    dbsnp="${index_path}/dbSnp.vcf.gz"
else 
    echo "Missing dbSNP File: ${index_path}/dbSnp.vcf.gz"
    usage
fi

if [[ -a "${index_path}/cosmic.vcf.gz" ]]
then
    cosmic=${index_path}/cosmic.vcf.gz
else 
    echo "Missing InDel File: ${index_path}/cosmic.vcf.gz"
    usage
fi
baseDir="`dirname \"$0\"`"

if [ $algo == 'strelka2' ]
  then
    module load strelka/2.8.3 manta/1.2.0 snpeff/4.3q vcftools/0.1.14
    mkdir manta strelka
    configManta.py --normalBam ${normal} --tumorBam ${tumor} --referenceFasta ${reffa} --runDir manta
    manta/runWorkflow.py -m local -j 8
    configureStrelkaSomaticWorkflow.py --normalBam ${mnormal} --tumorBam ${mtumor} --referenceFasta ${reffa} --targeted --indelCandidates manta/results/variants/candidateSmallIndels.vcf.gz --runDir strelka
    strelka/runWorkflow.py -m local -j 8
    vcf-concat strelka/results/variants/*.vcf.gz | vcf-annotate -n --fill-type -n |vcf-sort |java -jar $SNPEFF_HOME/SnpSift.jar filter "((FILTER = 'PASS') & (GEN[*].DP >= 10))" | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' |bgzip > ${pair_id}.strelka.vcf.gz
fi

if [ $algo == 'virmid' ]
  then 
    module load snpeff/4.3q virmid/1.2 vcftools/0.1.14
    virmid -R ${reffa} -D ${tumor} -N ${normal} -s ${cosmic} -t $SLURM_CPUS_ON_NODE -M 2000 -c1 10 -c2 10
    perl $baseDir/addgt_virmid.pl ${tumor}.virmid.som.passed.vcf
    perl $baseDir/addgt_virmid.pl ${tumor}.virmid.loh.passed.vcf
    vcf-concat *gt.vcf | vcf-sort | vcf-annotate -n --fill-type -n | java -jar $SNPEFF_HOME/SnpSift.jar filter '((NDP >= 10) & (DDP >= 10))' | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' | bgzip > ${pair_id}.virmid.vcf.gz
fi

if [ $algo == 'speedseq' ]
  then 
    module load snpeff/4.3q speedseq/20160506 vcftools/0.1.14
    speedseq somatic -q 10 -t $SLURM_CPUS_ON_NODE -o sssom ${reffa} ${normal} ${tumor}
    vcf-annotate -H -n --fill-type sssom.vcf.gz | java -jar $SNPEFF_HOME/SnpSift.jar filter '((QUAL >= 10) & (GEN[*].DP >= 10))' | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' |bgzip > ${pair_id}.sssom.vcf.gz
fi

if [ $algo == 'mutect2' ]
then
  module load parallel gatk/3.7 snpeff/4.3q vcftools/0.1.14
  cut -f 1 ${index_path}/genomefile.5M.txt | parallel --delay 2 -j 10 "java -Xmx20g -jar \$GATK_JAR -R ${reffa} -D ${dbsnp} -T MuTect2 -stand_call_conf 30 -stand_emit_conf 10.0 -A FisherStrand -A QualByDepth -A VariantType -A DepthPerAlleleBySample -A HaplotypeScore -A AlleleBalance -I:tumor ${tumor} -I:normal ${normal} --cosmic ${cosmic} -o ${tid}.{}.mutect.vcf -L {}"
  vcf-concat ${tid}*.vcf | vcf-sort | vcf-annotate -n --fill-type | java -jar \$SNPEFF_HOME/SnpSift.jar filter -p '((FS <= 60) & GEN[*].DP >= 10)' | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' |bgzip > ${pair_id}.pmutect.vcf.gz
fi

if [ $algo == 'varscan' ]
then
  module load snpeff/4.3q samtools/1.6 VarScan/2.4.2 speedseq/20160506 vcftools/0.1.14
  sambamba mpileup -t $SLURM_CPUS_ON_NODE ${tumor} --samtools "-C 50 -f ${reffa}"  > t.mpileup
  sambamba mpileup -t $SLURM_CPUS_ON_NODE ${normal} --samtools "-C 50 -f ${reffa}"  > n.mpileup
  VarScan somatic n.mpileup t.mpileup vscan --output-vcf 1
  VarScan copynumber n.mpileup t.mpileup vscancnv 
  vcf-concat vscan*.vcf | vcf-sort | vcf-annotate -n --fill-type -n | java -jar $SNPEFF_HOME/SnpSift.jar filter '((exists SOMATIC) & (GEN[*].DP >= 10))' | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' | bgzip > ${tid}_${nid}.varscan.vcf.gz
fi

if [ $algo == 'shimmer' ]
then
    module load snpeff/4.3q shimmer/0.1.1 vcftools/0.1.14
    shimmer.pl --minqual 25 --ref ${reffa} ${normal} ${tumor} --outdir shimmer 2> shimmer.err
    perl /project/PHG/PHG_Clinical/clinseq_workflows/scripts/add_readct_shimmer.pl
    vcf-annotate -n --fill-type shimmer/somatic_diffs.readct.vcf | java -jar $SNPEFF_HOME/SnpSift.jar filter '(GEN[*].DP >= 10)' | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' | bgzip > ${pair_id}.shimmer.vcf.gz
fi

if [ $algo == 'lancet' ]
then
    module load snpeff/4.3q lancet vcftools/0.1.14
    lancet --tumor ${tumor} --normal ${normal} --ref $reffa --bed $target_panel --num-threads 16 > out.vcf
    vcf-concat out.vcf | vcf-sort | vcf-annotate -n --fill-type -n | perl -pe 's/TUMOR/${tid}/' | perl -pe 's/NORMAL/${nid}/g' |bedtools intersect -header -a stdin -b $target_panel |bgzip > ${tid}_${nid}.lancet.vcf.gz
