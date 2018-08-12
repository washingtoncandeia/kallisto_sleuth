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

# Transcriptoma de referência para o index kallisto
REF=refs/GRCh38.cdna.fa
IDX=refs/GRCh38.cdna.fa.idx

# 1. Download do transcriptoma a partir do ftp Ensembl:
curl ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz | gunzip -c > $REF


# 2. Construção do index kallisto
kallisto index -i $IDX  $REF
