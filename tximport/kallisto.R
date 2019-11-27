## Kallisto
# Data: 26/11/2019

library(readr)
library(tximport)
library(biomaRt)
library(DESeq2)
library(readr)
library(dplyr)
library(rhdf5)
library(RColorBrewer)


## Parte 1 - Preparação de dados das amostras de kallisto.
# Caminho dos arquivos (fele path)
dir <- './results'
list.files(dir)


# Nomes de populações
ZIKA <- 'ZIKA'
CHIKV <- 'CHIIKV'
CHIKV_REC <- 'CHIKV_REC'
GBS <- 'GBS'
GBS_REC <- 'GBS_REC'
CONTROL <- 'CONTROL'

# Vetor com nomes populações para fazer coluna pop:
# Obs.: Ao fazer o vetor, melhor deixar organizado combinando amostra + nomearquivo
pop <- c(rep(GBS, 8), rep(CONTROL, 16), rep(ZIKA, 8), rep(GBS, 16), 
         rep(GBS_REC, 8), rep(CHIKV_REC, 8), rep(CONTROL, 8), rep(GBS, 8), rep(CONTROL, 8),
         rep(CHIKV, 8), rep(ZIKA, 8), rep(GBS, 8), rep(GBS_REC, 24), 
         rep(CHIKV, 16), rep(CONTROL, 16), rep(CHIKV_REC, 8), 
         rep(ZIKA, 8), rep(GBS_REC, 16), rep(GBS, 8), 
         rep(CHIKV, 8), rep(CHIKV_REC, 32), rep(ZIKA, 8),
         rep(GBS_REC, 24))

pop          # ZIKA, CHIKV, CHIKV_REC, GBS, GBS_REC, CONTROL
length(pop)  # 280

# Nome de centro de pesquisa para fazer coluna center
center <- rep('IMT', 280)
center
length(center) # 280 == pop

# Nomes de amostras para fazer coluna samples
amostras <- list.files(dir)
amostras
mode(amostras)
length(amostras)  # 280 == pop



## Parte 2 - Unir cada vetor formando colunas de um data frame:
data_fr <- data.frame(pop = pop, 
                      center = center, 
                      samples = amostras)
head(data_fr, 6)
str(data_fr)


# Observar cada coluna:
data_fr$pop
data_fr$center
data_fr$samples


# Salvar a tabela no formato .txt (tsv)
write.table(data_fr, 'amostras.txt', sep = '\t')


# Criar um vetor nomeado apontando os arquivos de quantificação.
# Estes arquivos têm seus nomes anotados em uma tabela (samples.txt).
samples <- read.table('amostras.txt', header = TRUE, row.names = 1)
head(samples)
samples$sample
mode(samples)


# Obtendo cada arquivo de replicata das amostras usadas em kallisto:
files <- file.path(dir, samples$sample, 'abundance.h5')
files


# Usando biomaRt
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                         dataset = "hsapiens_gene_ensembl", 
                         host="www.ensembl.org")

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", 
                                     "ensembl_gene_id"), 
                      mart = mart)

t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id)

head(t2g, 20)

# Usando tximport para kallisto
# Estimativa de contagens a partir de kallisto:
txi.kallisto <- tximport(files, 
                         type = 'kallisto',
                         tx2gene = t2g,
                         txOut = TRUE)





















