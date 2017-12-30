##--------------------------------------------
# Sleuth - RNA-seq Virus Zika
# Washington Candeia de Araujo
# Instituto de Medicina Tropical - IMD - UFRN
# Script Likelihood Ratio Test
##--------------------------------------------

library("sleuth")
library("cowplot")
library("dplyr")
library("ggplot2")
library("jsonlite")

## 1-Redirecionando para os diretorios de resultados
base_dir <- "/scratch/global/wcdaraujo/sleuth_z"
sample_id <- dir(file.path(base_dir,"results"))
sample_id
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results", id))
# Observando os caminhos
print(kal_dirs)

## 2-Criando tabela
s2c <- read.table(file.path(base_dir, "metadata", "metadata_1.txt"), header = TRUE, stringsAsFactors=FALSE)
head(s2c)

# Aqui, usar dplyr para descricao de relacoes entre as amostras
s2c <- dplyr::select(s2c, sample = samples, condition)
print(s2c)
# Construir uma coluna na tabela de metadados com os caminhos
s2c <- mutate(s2c, path = kal_dirs)
# Vamos ver os paths de cada arquivo
print(s2c)

## 3-Criando objeto: carregar o objeto kallisto
so <- sleuth_prep(s2c, ~condition, read_bootstrap_tpm = T, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
# Observar o modelo estatistico
models(so)

# Live
sleuth_live(so)


