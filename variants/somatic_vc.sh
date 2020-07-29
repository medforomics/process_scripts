#!/bin/bash
#run_somatic_caller.sh

usage(){
    echo "-h --Help documentation for run_somatic_caller.sh"
    echo "-a --Somatic Workflow Method: strelka2, virmid, speedseq, mutect2, varscan, shimmer, lancet"
    echo "-r  --Path to Reference Genome with the file genome.fa"
    echo "-n --Normal"
    echo "-t --Tumor"
    echo "-p -- PairID"
    echo "-x --NormalID"
    echo "-y --TumorID"
    echo "-i --NormalBAM used for Mantra in the case of UMI consensus"
    echo "-j --TumorBAM used for Mantra in the case of UMI consensus"
    echo "-b --TargetBed"
    echo "Example: bash somatic_vc.sh -a strelka2 -p subj1 -y ORD1_N_panel1385 -y ORD1_T_panel138 -n ORD1_N_panel1385.final.bam -t ORD1_T_panel1385.final.bam"
    exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :n:t:p:r:x:y:i:j:q:a:b:h opt
do
    case $opt in 
	r) index_path=$OPTARG;;
	x) tid=$OPTARG;;
	y) nid=$OPTARG;;
	p) pair_id=$OPTARG;;
	n) normal=$OPTARG;;
	t) tumor=$OPTARG;;
	a) algo=$OPTARG;;
	b) tbed=$OPTARG;; 
	q) pon==$OPTARG;; 
	h) usage;;
    esac
done

source /etc/profile.d/modules.sh
module load htslib/gcc/1.8 samtools/gcc/1.8 snpeff/4.3q vcftools/0.1.14
export PATH=/project/shared/bicf_workflow_ref/seqprg/bin:$PATH

shift $(($OPTIND -1))

#Check for mandatory options
if [[ -z $normal ]] || [[ -z $tumor ]] || [[ -z $algo ]]; then
    echo $normal $tumor $algo
    usage
fi 
NPROC=$SLURM_CPUS_ON_NODE
if [[ -z $NPROC ]] 
then 
    NPROC=`nproc`
fi

ponopt='';
if [[ -f $pon ]]
then
    ponopt="--pon $pon"
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

cat ${reffa}.fai |cut -f 1 |grep -v decoy |grep -v 'HLA' |grep -v alt |grep -v 'chrUn' |grep -v 'random' > intervals.txt
if [[ -n $tbed ]]
then
    awk '{print $1":"$2"-"$3}' $tbed > intervals.txt
fi
if [[ -z $tid ]]
then
    tid=`samtools view -H ${tumor} |grep '^@RG' |perl -pi -e 's/\t/\n/g' |grep ID |cut -f 2 -d ':'`
fi
if [[ -z $nid ]]
then
    nid=`samtools view -H ${normal} |grep '^@RG' |perl -pi -e 's/\t/\n/g' |grep ID |cut -f 2 -d ':'`
fi


if [ $algo == 'strelka2' ]
then
    module load strelka/2.9.10 manta/1.3.1 
    opt=''
    if [[ -n $tbed ]]
    then
	if [[ -f $tbed ]]
	then
	    opt="--callRegions ${tbed}.gz"
	else
	    cp $tbed targetpanel.bed
	    bgzip targetpanel.bed
	    tabix targetpanel.bed.gz
	    opt="--callRegions targetpanel.bed.gz"
	fi
    fi
    mkdir manta     
    configManta.py --normalBam ${normal} --tumorBam ${tumor} --referenceFasta ${reffa} $opt --runDir manta
    manta/runWorkflow.py -m local -j 8
    mantaopt=''
    if [[ -f manta/results/variants/candidateSmallIndels.vcf.gz ]]
    then
	mantaopt="--indelCandidates manta/results/variants/candidateSmallIndels.vcf.gz"
    fi
    configureStrelkaSomaticWorkflow.py --normalBam ${normal} --tumorBam ${tumor} --referenceFasta ${reffa} --targeted --runDir strelka $mantaopt
    strelka/runWorkflow.py -m local -j 8
    vcf-concat strelka/results/variants/*.vcf.gz | vcf-annotate -n --fill-type -n |vcf-sort |java -jar $SNPEFF_HOME/SnpSift.jar filter "(GEN[*].DP >= 10)" | perl -pe "s/TUMOR/${tid}/g" | perl -pe "s/NORMAL/${nid}/g" |bgzip > ${pair_id}.strelka2.vcf.gz
elif [ $algo == 'virmid' ]
then 
    module load virmid/1.2
    virmid -R ${reffa} -D ${tumor} -N ${normal} -s ${cosmic} -t $NPROC -M 2000 -c1 10 -c2 10
    perl $baseDir/addgt_virmid.pl ${tumor}.virmid.som.passed.vcf
    perl $baseDir/addgt_virmid.pl ${tumor}.virmid.loh.passed.vcf
    module rm java/oracle/jdk1.7.0_51
    module load snpeff/4.3q
    vcf-concat *gt.vcf | vcf-sort | vcf-annotate -n --fill-type -n | java -jar $SNPEFF_HOME/SnpSift.jar filter '((NDP >= 10) & (DDP >= 10))' | perl -pe "s/TUMOR/${tid}/g" | perl -pe "s/NORMAL/${nid}/g" | bgzip > ${pair_id}.virmid.vcf.gz
elif [ $algo == 'mutect' ]
then
    gatk4_dbsnp=${index_path}/clinseq_prj/dbSnp.gatk4.vcf.gz
    module load gatk/4.1.4.0 parallel/20150122
    threads=`expr $NPROC / 2`
    cut -f 1 intervals.txt | parallel --delay 1 --jobs $threads "gatk --java-options \"-Xmx20g\" Mutect2 $ponopt --independent-mates -RF AllowAllReadsReadFilter -R ${reffa} -I ${tumor} -tumor ${tid} -I ${normal} -normal ${nid} --output ${tid}.mutect.{}.vcf -L {}"
    vcf-concat ${tid}.mutect.*vcf | vcf-sort | vcf-annotate -n --fill-type | java -jar $SNPEFF_HOME/SnpSift.jar filter -p '(GEN[*].DP >= 10)' | bgzip > ${pair_id}.mutect.vcf.gz
elif [ $algo == 'varscan' ]
then
    module load bcftools/gcc/1.8 VarScan/2.4.2
    module rm java/oracle/jdk1.7.0_51
    module load snpeff/4.3q 
    samtools mpileup -C 50 -f ${reffa} $tumor > t.mpileup
    samtools mpileup -C 50 -f ${reffa} $normal > n.mpileup
    VarScan somatic n.mpileup t.mpileup vscan --output-vcf 1
    VarScan copynumber n.mpileup t.mpileup vscancnv
    vcf-concat vscan*.vcf | vcf-sort | vcf-annotate -n --fill-type -n | java -jar $SNPEFF_HOME/SnpSift.jar filter '((exists SOMATIC) & (GEN[*].DP >= 10))' | perl -pe "s/TUMOR/${tid}/" | perl -pe "s/NORMAL/${nid}/g" | bgzip >  ${pair_id}.varscan.vcf.gz
elif [ $algo == 'shimmer' ]
then
    module load R/3.6.1-gccmkl
    shimmer.pl --minqual 25 --ref ${reffa} ${normal} ${tumor} --outdir shimmer 2> shimmer.err
    perl $baseDir/add_readct_shimmer.pl
    module rm java/oracle/jdk1.7.0_51
    module load snpeff/4.3q
    vcf-annotate -n --fill-type shimmer/somatic_diffs.readct.vcf | java -jar $SNPEFF_HOME/SnpSift.jar filter '(GEN[*].DP >= 10)' | perl -pe "s/TUMOR/${tid}/" | perl -pe "s/NORMAL/${nid}/g" | bgzip > ${pair_id}.shimmer.vcf.gz
fi


