#!/bin/bash
#SBATCH --job-name=kallisto 
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=17-0:0


# Protocolo: Kallisto RNA-seq
# Autor: Washington Candeia de Araujo
# Data: 10 Ago 2017

# Criar diretório para armazenar resultados
mkdir -p results info 

# Encerrar script se houver algum erro.
set -euo pipefail


echo "Análise feita por `whoami` em `date`" 

# Transcriptoma de referência para o index kallisto
REF=refs/Homo_sapiens.GRCh38.rna.fa.gz 

# O nome do index gerado por kallisto
IDX=refs/hsGRCh38_kallisto.idx

# Diretório de armazenamento dos resultados:
OUTDIR=results/

echo "INICIANDO ANÁLISE COM KALLISTO"

# Kallisto Quant
for SAMPLE in {1..48} 
do 
  for REPLICATE in {1..8}
  do 
     R1=new_selma/samples/${SAMPLE}_lane${REPLICATE}_20161227000_L00${REPLICATE}_R1_001.fastq.gz;
     R2=new_selma/samples/${SAMPLE}_lane${REPLICATE}_20161227000_L00${REPLICATE}_R2_001.fastq.gz; 

     echo "Kallisto quantificação: ";
     echo "Amostras (paired end): $R1 + $R2";
            
     kallisto quant -i $IDX -t 32 -b 100 $R1 $R2 -o $OUTDIR${SAMPLE}_rep${REPLICATE}; 
  done
done; 




