#!/bin/bash

processName=$1
outputDir=./
inputDir=./
scriptsDir=/home/sbsuser/Scritps/benchmarking/

grep "Elapsed time" ${inputDir}${processName}.stdout | python ${scriptsDir}extractIsisSteps.py > $outputDir$processName

python ${scriptsDir}collectlTimeDateConverter.py "$processName" "$outputDir" "$inputDir"
cp ${processName}.tsv collectl.tsv

echo $outputDir

echo $filename
export GNUTERM=dumb
/usr/bin/gnuplot <<\EOF

reset
clear
workflow = ""
hostname = ""
file="collectl.tsv"
set terminal png size 1500,500
set output "CPU".hostname.workflow.".png"

set xdata time
set timefmt "%s"
set format x "%H:%M"
set xlabel "time"
set ylabel "Percent CPU"
#set yrange [0:110]
set title "CPU usage"." ".hostname." ".workflow
set style data line
set key outside
#set xtic (1000)

plot file using (column("Time")):(column("[CPU]Wait%")):xtic(60) title "[CPU]Wait%", "" using (column("Time")):(column("[CPU]Sys%")) title "[CPU]Sys%", "" using (column("Time")):(column("[CPU]Nice%")) title "[CPU]Nice%", "" using (column("Time")):(column("[CPU]User%")) title "[CPU]User%"

set output "memory".hostname.workflow.".png"
set xdata time
set timefmt "%s"
set format x "%H:%M"
set xlabel "time"
set ylabel "GB"
set title "Memory usage"." ".hostname." ".workflow
set style data line

plot file using (column("Time")):(column("[MEM]Commit")) title "[MEM]Commit", "" using (column("Time")):(column("[MEM]Cached")) title "[MEM]Cashed"

set output "IO".hostname.workflow.".png"
set xdata time
set timefmt "%s"
set format x "%H:%M"
set xlabel "time"
set ylabel "MB/s
set title "Disk IO"." ".hostname." ".workflow
set style data line
plot file using (column("Time")):(column("[DSK]ReadKBTot")) title "[DSK]ReadKBTot%", "" using (column("Time")):(column("[DSK]WriteKBTot")) title "[DSK]WriteKBTot%"
set output
#
EOF
rm collectl.tsv
