#!/bin/bash

outDir='../Out/'
for i in `seq 0 2 100`
do  
  echo $i
  outFile=$outDir"gl_"$i"_.RData"
  echo $outFile
  /Library/Frameworks/R.framework/Resources/Rscript glassoScript.R $1 $2 $3
done