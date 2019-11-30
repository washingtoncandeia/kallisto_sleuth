## Análise em Nível de Transcrito Kallisto
## Utilizando Triplicatas
# Data: 30/11/2019

library(tximport)
library(biomaRt)
library(DESeq2)
library(readr)
library(dplyr)
library(rhdf5)
library(RColorBrewer)
library(dplyr)

## Parte 1 - Preparação de dados das amostras de kallisto.
# Caminho dos arquivos (fele path)
dir <- './results'
list.files(dir)


# Nomes de populações
ZIKA <- 'ZIKA'
CHIKV <- 'CHIKV'
CHIKV_REC <- 'CHIKV_REC'
GBS <- 'GBS'
GBS_REC <- 'GBS_REC'
CONTROL <- 'CONTROL'

# Vetor com nomes populações para fazer coluna pop:
# Obs.: Ao fazer o vetor, melhor deixar organizado combinando amostra + nomearquivo
pop <- c(rep(GBS, 3), rep(CONTROL, 6), rep(ZIKA, 3), rep(GBS, 6), 
         rep(GBS_REC, 3), rep(CHIKV_REC, 3), rep(CONTROL, 3), rep(GBS, 3), 
         rep(CONTROL, 3),rep(CHIKV, 3), rep(ZIKA, 3), rep(GBS, 3), 
         rep(GBS_REC, 9), rep(CHIKV, 6), rep(CONTROL, 6), rep(CHIKV_REC, 3), 
         rep(ZIKA, 3), rep(GBS_REC, 6), rep(GBS, 3), rep(CHIKV, 3), 
         rep(CHIKV_REC, 12), rep(ZIKA, 3), rep(GBS_REC, 9))

pop              # ZIKA, CHIKV, CHIKV_REC, GBS, GBS_REC, CONTROL
length(pop)      # 105

# Nome de centro de pesquisa para fazer coluna center
center <- rep('IMT', 105)
center
length(center)   # 105 == pop

# Nomes de amostras para fazer coluna samples
amostras <- list.files(dir)
amostras
mode(amostras)
length(amostras)  # 105 == pop



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

## -------------------- Eliminando Linhas do Data Frame por Nomes  -------------------- ##
## Quando necessário o uso de um data frame menor, com apenas algumas variáveis.

# Usando dplyr, função filter e negando com regex (função grepl)

data_fr <- data_fr %>% 
  filter(!grepl('GBS', pop))

data_fr <- data_fr %>% 
  filter(!grepl('GBS_REC', pop))

data_fr <- data_fr %>% 
  filter(!grepl('CHIKV', pop))

data_fr <- data_fr %>% 
  filter(!grepl('CHIKV_REC', pop))

data_fr

## ------------------------------------------------------------------------------------ ##


# Salvar a tabela no formato .txt (tsv)
write.table(data_fr, 'amostras_zikvCtl.txt', sep = '\t')


# Criar um vetor nomeado apontando os arquivos de quantificação.
# Estes arquivos têm seus nomes anotados em uma tabela (samples.txt).
samples <- read.table('amostras_zikvCtl.txt', header = TRUE, row.names = 1)
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



# Salvamento de um objeto R para uso posterior:
dir.create(path = "./count_estimates/")
save(txi.kallisto, file = "./count_estimates/txi_count_Trip_estimates.Rdata")


# Observar a lista gerada
head(txi.kallisto)

## Passo 4 - Usando DESeq2 para anaĺise de Expressão Diferencial (DE)
samples_information <- data_fr
head(samples_information)

samples_information

samples_information$pop <- relevel(x = samples_information$pop, ref = CONTROL)



## Wald test para DE:
# Quais genes estão diferencialmente expressos entre cada tratamento e controles.
## Define reference and alternative caste:
reference_caste <- CONTROL
alternative_caste <- ZIKA
castes <- c(reference_caste, alternative_caste)



levels(samples_information)

## Set ZIKV as reference for caste:
samples_information$pop <- relevel(x = samples_information$pop,
                                     ref = reference_caste)


deseq_txi <- DESeqDataSetFromTximport(txi = txi.kallisto,
                                      colData = samples_information,
                                      design = ~pop + pop)
deseq_txi


## Perform a pairwise comparison between treatments:
deseq_object  <- DESeq(deseq_txi,
                       test = "Wald",
                       betaPrior = FALSE)

## Extract information related to treatments:
treatments <- c(ZIKA, CONTROL)


