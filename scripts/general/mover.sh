#!/bin/bash
# 12 de agosto de 2018

DIR=/data/home/waraujo/Selma
COPY=/data/home/waraujo/samples/

# Utilizando for loop para mover todos os arquivos.

for SAMPLE in {1..48}
do
  for REPLICATE in {1..8}
  do
    R1=${DIR}/Sample_${SAMPLE}_lane${REPLICATE}_20161227000/${SAMPLE}_lane${REPLICATE}_20161227000_L00${REPLICATE}_R1_001.fastq.gz;
    R2=${DIR}/Sample_${SAMPLE}_lane${REPLICATE}_20161227000/${SAMPLE}_lane${REPLICATE}_20161227000_L00${REPLICATE}_R2_001.fastq.gz;
    
    echo "Amostras R1: $R1";
    echo "Amostras R2: $R2";
    
    mv ${R1} ${SAMPLES};
    mv ${R2} ${SAMPLES};
  done;
done;
