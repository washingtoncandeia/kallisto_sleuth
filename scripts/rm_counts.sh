#!/bin/bash

# Remover arquivos iniciados com _

DIR=/home/wash/kallisto/results   

cd $DIR
# Removendo
for SAMPLE in {1..48}_;
do
  for REPLICATE in rep{1..8};
  do
    rm ${DIR}/${SAMPLE}${REPLICATE}/_${REPLICATE}.counts.tsv
  done;
done;
