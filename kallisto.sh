#!/bin/bash
##----------------------------------------------
# Protocolo: Kallisto RNA-seq
# Autor: Washington Candeia de Araujo
# Data modificação: 12 de agosto de 2018
##-----------------------------------------------

# Criar diretório para armazenar resultados
mkdir -p results

# Encerrar script se houver algum erro.
set -euo pipefail

# Diretório amostras
DIR=/data/home/waraujo/samples

# Transcriptoma de referência para o index kallisto
REF=/data/home/waraujo/kallisto/refs/GRCh38.cdna.fa

# Index kallisto
IDX=/data/home/waraujo/kallisto/refs/GRCh38.cdna.fa.idx

# Diretório de armazenamento dos resultados:
OUTDIR=./results/


# 1. Kallisto Quant
for SAMPLE in {1..48} 
do 
  for REPLICATE in {1..8}
  do 
     R1=${DIR}/${SAMPLE}_lane${REPLICATE}_20161227000_L00${REPLICATE}_R1_001.fastq.gz;
     R2=${DIR}/${SAMPLE}_lane${REPLICATE}_20161227000_L00${REPLICATE}_R2_001.fastq.gz; 

     echo "Kallisto quantificação: ";
     echo "Amostra: $R1";
     echo "Amostra: $R2";      
     
     kallisto quant -i $IDX -t 64 -b 100 $R1 $R2 -o $OUTDIR${SAMPLE}_rep${REPLICATE}; 
  done;
done; 

