#!/bin/bash
#trimgalore.sh

pair_id=$1
index_path=$2
fq1=$3
fq2=$4

module add python/2.7.x-anaconda star/2.5.2b
STAR-Fusion --genome_lib_dir ${index_path}/CTAT_lib/ --left_fq ${fq1} --right_fq ${fq2} --output_dir star_fusion &> star_fusion.err
mv star_fusion/star-fusion.fusion_candidates.final.abridged ${pair_id}.starfusion.txt
