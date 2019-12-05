#!/bin/bash

# Data: 03/08/2017

## Renomear os arquivos do RNA-seq
## retirando o [S.+] ou _[S]\d+_

for nome in `ls */*.fastq.gz`; do
    novo=`echo $nome | sed -r 's/(.+)_S.+_(L.+)/\1_\2/'`;
    echo $novo;
    mv $nome $novo
done;
