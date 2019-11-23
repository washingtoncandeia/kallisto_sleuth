## Kallisto
# Data: 19/11/2019

library(readr)
library(tximport)

# Caminho dos arquivos (fele path)
dir <- './results'
list.files(dir)

# Criar um vetor nomeado apontando os arquivos de quantificação.
# Estes arquivos têm seus nomes anotados em uma tabela (samples.txt).
samples <- read.table('samples.txt', header = TRUE)
head(samples)
samples$sample


files <- file.path(dir, samples$sample, 'abundance.h5')
files

txi.kallisto <- tximport(files, type = 'kallisto', txOut = TRUE)
