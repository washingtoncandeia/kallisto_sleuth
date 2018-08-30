##--------------------------------------------
# Sleuth - RNA-seq Virus Zika
# Washington Candeia de Araujo
# Instituto de Medicina Tropical - IMD - UFRN
# Script lrt
# https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html
# Data: 30/08/2018
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
# metadata_2.txt (GBS-rec x ZIKV x control).
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

## 3-Criando objeto: carregar o objeto kallisto.
#so <- sleuth_prep(s2c, ~condition, read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE, target_mapping = t2g)
#-- extra_bootstrap_summary: se TRUE, computa estatística extra para contagens estimadas. Necessária para plot_bootstrap.

#so <- sleuth_prep(s2c, ~condition, read_bootstrap_tpm = TRUE, max_bootstrap = 30, target_mapping = t2g)
#-- max_bootstrap: o número máximo de valores de bootstrap para ler para cada transcrito. Baixar diminui a acurácia da estimativa
#                  do ruído inferencial.
so <- sleuth_prep(s2c, ~condition, target_mapping = t2g)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

# Observar o modelo estatistico
models(so)

# 4. Live (opcional).
sleuth_live(so)

# 5. Gerar resultados:
sleuth_table <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
#sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

# 5.1. Filtrar por valores significantes.
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

# 5.2. Obter genes significantes.
head(sleuth_significant, 20)
# Esta tablea mostra os 20 genes mais significantes, 
# com um q-value <= 0.05 (Benjamini-Hochberg multiple testing corrected).

# 5.3. Escrever uma tabela com resultados de sleuth.
write.table(as.data.frame(head(sleuth_significant, 400)), sep = "\t", "zika_significant_1.txt", row.names = F)

# 5.4. Bootstrap.
#plot_bootstrap(so, "ENST00000263734", units = "est_counts", color_by = "condition")

# 5.5. Plot PCA.
plot_pca(so, text_labels = TRUE, color_by = 'sample')

# 5.6. Controle de qualidade.
plot_group_density(so, use_filtered = TRUE, units = "est_counts",
                   trans = "log", grouping = setdiff(colnames(so$sample_to_covariates), 
                   "sample"), offset = 1)
