#!/bin/bash

# Copiar diretórios de resultados das análises
# saídas de kallisto.
# Controles (10,11,19,20,33,34)

DIR=/scratch/global/wcdaraujo/results

# Utilizando for loop para copiar recursivamente.
# /scratch/global/wcdaraujo/results
for SAMPLE in 12 24 36 48;
do
  for REPLICATE in _rep{1..3};
  do
    cp -r ${DIR}/${SAMPLE}${REPLICATE} . 
  done;
done;

