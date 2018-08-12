#!/bin/bash

# Washington Candeia de Araújo
# Data: 20/08/2017
# Copiar todos os arquivos abundance.tsv, para o diretório ./copy_abundance
# renomeando as cópias.

DIR=/scratch/global/wcdaraujo/results/ 
COPY=./copy_abundance/


# Criando diretórios para todos os arquivos abundance.tsv
for sample in {1..48}_;
do
  for REPLICATE in rep{1..8};
  do
     mkdir ${COPY}${sample}${REPLICATE}
  done;
done;

# Copiando todos os arquivos abundance.tsv para cada diretório.
# Renomeando-os conforme a amostra.
for sample in {1..48}_;
do
  for REPLICATE in rep{1..8};
  do
     cp ${DIR}${sample}${REPLICATE}/abundance.tsv ${COPY}${sample}${REPLICATE}/${sample}${REPLICATE}.counts.tsv;
  done;
done;

