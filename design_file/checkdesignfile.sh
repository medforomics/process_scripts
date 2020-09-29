#!/bin/bash
#check_inputfiles.sh

baseDir="`dirname \"$0\"`"

rpair=$1
design=$2

perl -p -e 's/\\r\\n*/\\n/g' $design > design.fix.txt
perl $baseDir/check_designfile.pl ${rpair} design.fix.txt
