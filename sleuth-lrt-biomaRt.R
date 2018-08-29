##--------------------------------------------
# Sleuth - RNA-seq Virus Zika
# Washington Candeia de Araujo
# Versão GNULinux (29/08/2018)
# Instituto de Medicina Tropical - IMD - UFRN
# Script lrt
# https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html
##--------------------------------------------

library("sleuth")
library("cowplot")
library("dplyr")
library("ggplot2")
library("jsonlite")
library("biomaRt")

## 1-Redirecionando para os diretorios de resultados
base_dir <- "~/kallisto"
sample_id <- dir(file.path(base_dir,"results"))
sample_id
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results", id))
# Observando os caminhos
print(kal_dirs)

## 2-Criando tabela
# metadata_2 (GBS-rec x ZIKV x control)
s2c <- read.table(file.path(base_dir, "metadata", "metadata_2.txt"), header = TRUE, stringsAsFactors = FALSE)
head(s2c)

# Aqui, usar dplyr para descricao de relacoes entre as amostras
s2c <- dplyr::select(s2c, sample = samples, condition)
print(s2c)
# Construir uma coluna na tabela de metadados com os caminhos
s2c <- mutate(s2c, path = kal_dirs)
# Vamos ver os paths de cada arquivo
print(s2c)

## 2.1 - Usando biomaRt aqui: https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)


## 3-Criando objeto: carregar o objeto kallisto
so <- sleuth_prep(s2c, ~condition, target_mapping = t2g)
#so <- sleuth_prep(s2c, ~condition, read_bootstrap_tpm = T, extra_bootstrap_summary = TRUE)
#so <- sleuth_prep(s2c, ~condition, read_bootstrap_tpm = T, max_bootstrap = 30)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

# Observar o modelo estatistico
models(so)

sleuth_live(so)

# Gerar resultados:
results_table <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
#sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

# Obter genes significantes:
head(sleuth_significant, 20)
# The table shown above displays the top 20 significant genes 
# with a (Benjamini-Hochberg multiple testing corrected) q-value <= 0.05

# Escrever uma tabela com resultados de sleuth
write.table(as.data.frame(head(sleuth_significant, 400)), sep = "\t", "zika_significant_1.txt", row.names = F)

#plot_pca(so, text_labels = TRUE, color_by = 'sample')
