#!/bin/bash
##----------------------------------------------
# Protocolo: Kallisto RNA-seq
# Autor: Washington Candeia de Araujo
# Data modificação: 13 de agosto de 2018
##-----------------------------------------------

# Encerrar sob qualquer erro
set -p euo pipefail

# Criar um diretório para armazenar resultados
mkdir -p results

# Criar um diretório para armazenar indice e trasncriptoma referência
mkdir -p refs

# Diretório amostras
DIR=/data/home/waraujo/samples

# Transcriptoma de referência para o index kallisto
REF=/data/home/waraujo/kallisto/refs/GRCh38.cdna.fa

# Index kallisto
IDX=/data/home/waraujo/kallisto/refs/GRCh38.cdna.fa.idx

# Diretório de armazenamento dos resultados:
OUTDIR=./results/

# 1. Download do transcriptoma a partir do ftp Ensembl:
if [ ! -f $REF ]; then
    
    # A URL do transcriptoma
    REF_URL=ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz 
    
    echo "*** Baixar o transcriptoma de referência: $REF_URL"
    curl $REF_URL | gunzip -c > $REF
else
    echo "*** Transcriptoma de referência encontrado: $REF."
fi


# 2. Construção do index kallisto
kallisto index -i $IDX $REF

# 3. Kallisto Quant
for SAMPLE in {1..48} 
do 
  for REPLICATE in {1..8}
  do 
    R1=${DIR}/${SAMPLE}_lane${REPLICATE}_20161227000_L00${REPLICATE}_R1_001.fastq.gz;
    R2=${DIR}/${SAMPLE}_lane${REPLICATE}_20161227000_L00${REPLICATE}_R2_001.fastq.gz; 

    echo "Kallisto quantificação: ";
    echo "Amostra: $R1";
    echo "Amostra: $R2";      
     
    kallisto quant -i $IDX -t 16 -b 100 $R1 $R2 -o $OUTDIR${SAMPLE}_rep${REPLICATE}; 
  done;
done;
