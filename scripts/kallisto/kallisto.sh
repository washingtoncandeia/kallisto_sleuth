#!/bin/bash
##----------------------------------------------
# Protocolo: Kallisto RNA-seq
# Autor: Washington Candeia de Araujo
# Data modificação: 11 de agosto de 2018
# Modificações feitas no download ftp
# e na $REF usado
##-----------------------------------------------

# Criar diretório para armazenar resultados
mkdir -p results info refs

# Encerrar script se houver algum erro.
set -euo pipefail

GTF=refs/GRCh38.gtf
# Download GTF
if [ ! -f $GTF ]; then
 
    # Sitio GTF Ensembl:
    GTF_URL=ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz
 
    echo "*** Downloading the GTF from: $GTF"
    curl $GTF_URL | gunzip -c > $GTF
else
     echo "*** Found GTF file: $GTF"
fi

# Transcriptoma de referência para o index kallisto
REF=refs/GRCh38.cdna.fa
IDX=refs/GRCh38.cdna.fa.idx

# 1. Download do transcriptoma a partir do ftp Ensembl:
curl ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz | gunzip -c > $REF


# 2. Construção do index kallisto
kallisto index -i $IDX  $REF

# 3. Kallisto Quant
for SAMPLE in {1..48}; 
do 
  for REPLICATE in {1..8}
  do 
     R1=Project_Selma/${SAMPLE}_lane${REPLICATE}_20161227000_L00${REPLICATE}_R1_001.fastq.gz;
     R2=Project_Selma/${SAMPLE}_lane${REPLICATE}_20161227000_L00${REPLICATE}_R2_001.fastq.gz; 

     echo "Kallisto quantificação: ";
     echo "Amostras (paired end): $R1 + $R2";
      
     $OUTDIR=results/${SAMPLE}_lane${REPLICATE}
     kallisto quant -i $IDX -t 32 -b 100 $R1 $R2 -o $OUTDIR${SAMPLE}_rep${REPLICATE}; 
  done;
done; 




