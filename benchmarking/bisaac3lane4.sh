#!/bin/bash

system=`hostname`
workflow=$1
run=$2
processName=${system}${workflow}${run}
workingdir=./
analysisdir=${workingdir}runs/$processName/
command=/illumina/development/iSAAC/iSAAC-03.15.07.30/bin/isaac-align
validationFile=${analysisdir}Align/Projects/default/default/sorted.bam.md5
resultsDir=${workingdir}results/

echo "this script will run iSAAC benchmarking and monitor system using collectl"

mkdir $analysisdir

collectl -s cmd -f $analysisdir$processName &
CollectlPid=$!

echo collectl started

nohup /usr/bin/time -v $command -r /illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/IsaacIndex5/sorted-reference.xml -b RunInfo.xml --base-calls-format bcl-gz -m 122 --tiles s_4 -o ${analysisdir}/Align --buffer-bins on > ${workingdir}${analysisdir}${processName}.stdout

echo killing collectl

kill $CollectlPid

collectl -p $analysisdir${processName}* -P -f ${analysisdir}plot*

gzip -d ${analysisdir}plot*

mv ${analysisdir}plot* ${analysisdir}${processName}.dat

cp $validationFile ${resultsDir}${processName}.val


echo $processName >> ${resultsDir}runtimes
grep "Elapsed (wall clock) time" ${analysisdir}${processName}.stdout >> ${resultsDir}runtimes

exit
