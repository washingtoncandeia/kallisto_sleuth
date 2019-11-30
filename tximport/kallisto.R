## Análise em Nível de Transcrito Kallisto
## Utilizando Triplicatas
# Data: 30/11/2019

library(tximport)
library(biomaRt)
library(DESeq2)
library(readr)
library(dplyr)
library(rhdf5)
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
center <- rep('IMT-UFRN', 105)
center
length(center)   # 105 == pop


# Nomes de amostras analisadas para fazer coluna run
run <- list.files(dir)
run
mode(run)
length(run)  # 105 == pop


condition <- c(rep('gbs', 3), rep('control', 6), rep('zika', 3), rep('gbs', 6), 
               rep('gbs_rec', 3), rep('chikv_rec', 3), rep('control', 3), rep('gbs', 3), 
               rep('control', 3),rep('chikv', 3), rep('zika', 3), rep('gbs', 3), 
               rep('gbs_rec', 9), rep(CHIKV, 6), rep('control', 6), rep('chikv_rec', 3), 
               rep('zika', 3), rep('gbs_rec', 6), rep('gbs', 3), rep(CHIKV, 3), 
               rep('chikv_rec', 12), rep('zika', 3), rep('gbs_rec', 9))


# Replicatas
replicates <- c('rep01', 'rep02', 'rep03 ')


## Parte 2 - Unir cada vetor formando colunas de um data frame:
samples_info <- data.frame(pop = pop, 
                      center = center, 
                      run = run,
                      condition = condition,
                      replicate = rep(replicates, 35))   

# Obs.: Replicate: Não colocar 105 pois são de 3 em 3.
# Logo: 35 x 3 = 105. Repete-se a tríade 35 vezes, o que geram 105 replicatas.


head(samples_info, 10)
str(samples_info)
names(samples_info)


# Observar cada coluna:
samples_info$pop
samples_info$center
samples_info$run
samples_info$condition
samples_info$replicate
str(samples_info$pop)
str(samples_info$condition)
str(samples_info$run)
str(samples_info$replicate)


## -------------------- Eliminando Linhas do Data Frame por Nomes  -------------------- ##
## Quando necessário o uso de um data frame menor, com apenas algumas variáveis.

# Usando dplyr, função filter e negando com regex (função grepl)

samples_info <- samples_info %>% 
  filter(!grepl('GBS', pop))

samples_info <- samples_info %>% 
  filter(!grepl('GBS_REC', pop))

samples_info <- samples_info %>% 
  filter(!grepl('CHIKV', pop))

samples_info <- samples_info %>% 
  filter(!grepl('CHIKV_REC', pop))

samples_info

## ------------------------------------------------------------------------------------ ##


# Salvar a tabela no formato .txt (tsv)
write.table(samples_info, 'amostras_zikvCtl.txt', sep = '\t')


# Criar um vetor nomeado apontando os arquivos de quantificação.
# Estes arquivos têm seus nomes anotados em uma tabela (samples.txt).
samples <- read.table('amostras_zikvCtl.txt', header = TRUE, row.names = 1)
head(samples)
samples$condition 
mode(samples)

# Nomeando as linhas com nome de cada arquivo de amostra:
rownames(samples) <- samples$run

samples
samples[ , c('pop', 'center', 'run', 'condition', 'replicate')]


# Obtendo cada arquivo de replicata das amostras usadas em kallisto:
files <- file.path(dir, samples$run, 'abundance.h5')
files
names(files) <- samples$run
files

# Usando biomaRt
mart <- biomaRt::useMart(biomart = "ensembl", 
                         dataset = "hsapiens_gene_ensembl", 
                         host="www.ensembl.org")


t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", 
                                     "ensembl_gene_id"), 
                      mart = mart)


t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id)


head(t2g, 20)
## Obs.: A coluna de transcripts IDs não possui versão.
## Ao utilizar o tximport prestar atenção na opção ignoreTxVersion.


### Usando tximport para kallisto


# Estimativa de contagens a partir de kallisto,
# Usar ignoreTxVersion e ignoreAfterBar para que o data frame de IDs de transcritos
# e genes do Ensembl tenham ignorados as versões e barras |, respectivamente.
txi.kallisto <- tximport(files, 
                         type = 'kallisto',
                         tx2gene = t2g,
                         ignoreTxVersion = TRUE,
                         ignoreAfterBar = TRUE)

head(txi.kallisto)

# Salvamento de um objeto R para uso posterior:
dir.create(path = "./count_estimates/")
save(txi.kallisto, file = "./count_estimates/txi_count_Trip_estimates.Rdata")


# Observar a lista gerada
head(txi.kallisto)
names(txi.kallisto)
head(txi.kallisto$abundance)


#### DESeq2
## Design com formula simples:
dds.txi <- DESeqDataSetFromTximport(txi = txi.kallisto,
                                    colData = samples,
                                    design= ~ condition)

# Agora, o objeto dds.Txi pode ser usado como aquele dds nos
# passos subsequentes de DESeq2.


## Design com formula ~replicate + condition
# dds.txi <- DESeqDataSetFromTximport(txi = txi.kallisto,
#                                     colData = samples,
#                                     design= ~ replicate + condition)


## Observação de informarção após esse comando acima:
# https://support.bioconductor.org/p/62245/

# Michael Love responde:
# "Because you have other variables in the design other than Condition (Sex and ERCCMix), 
# you do not have 3 replicates of a unique combination of the design variables. 
# The count outlier behavior is explained in the vignette section "Approach to count outliers":

# "The results function automatically flags genes which contain a Cook’s distance above a cutoff 
# for samples which have 3 or more replicates."

# "With less than 3 replicates per unique combination, there is no filtering of potential count outliers."
