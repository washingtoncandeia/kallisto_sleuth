## Kallisto
# Data: 21/11/2019
# Criação de data frame amostras.txt do zero
# amostras.txt == samples.txt
##
library(readr)
library(tximport)

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

# Unir cada vetor formando colunas de uma tabela:
data_fr <- data.frame(pop = pop, center = center, samples = amostras)
head(data_fr, 6)
str(data_fr)

data_fr$pop
data_fr$center
data_fr$samples


write.table(data_fr, 'amostras.txt', sep = '\t')

# Criar um vetor nomeado apontando os arquivos de quantificação.
# Estes arquivos têm seus nomes anotados em uma tabela (samples.txt).
samples <- read.table('amostras.txt', header = T)
head(samples)
samples$sample
mode(samples)

files <- file.path(dir, samples$sample, 'abundance.h5')
files

txi.kallisto <- tximport(files, type = 'kallisto', txOut = TRUE)
