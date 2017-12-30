##--------------------------------------------
# Sleuth - RNA-seq Virus Zika
# Washington Candeia de Araujo
# Versão GNULinux (21/11/2017)
# Instituto de Medicina Tropical - IMD - UFRN
# Script Wald Test
# Supercomputador - Wald Test
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
s2c <- read.table(file.path(base_dir, "metadata", "metadata_6.txt"), header = TRUE, stringsAsFactors=FALSE)
head(s2c)

# Aqui, usar dplyr para descricao de relacoes entre as amostras
s2c <- dplyr::select(s2c, sample = samples, condition)
print(s2c)
# Construir uma coluna na tabela de metadados com os caminhos
s2c <- mutate(s2c, path = kal_dirs)
# Vamos ver os paths de cada arquivo
print(s2c)

## 3-Criando objeto: carregar o objeto kallisto
so <- sleuth_prep(s2c, ~condition, read_bootstrap_tpm = T, max_bootstrap = 30)
so <- sleuth_fit(so)
so <- sleuth_wt(so, which_beta = 'conditionctl')
# Observar o modelo estatistico
models(so)


##----------------------------------------------------------------
# Examinar os resultados do teste
sleuth_table <- sleuth_results(so, 'conditionctl', test_type = 'wt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
# Obter genes significantes:
head(sleuth_significant, 20)
# The table shown above displays the top 20 significant genes 
# with a (Benjamini-Hochberg multiple testing corrected) q-value <= 0.05
write.table(as.data.frame(head(sleuth_significant, 400)), sep = "\t", "significant_wt.txt", row.names = F)
##----------------------------------------------------------------

# Live
sleuth_live(so)

##----------------------------------------------------------------
# Tabela Nome dos Genes
# Examinar os resultados do teste, com resultados de isoformas:
sleuth_table <- sleuth_results(so, 'conditionctl', test_type = 'wt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
# Obter genes significantes:
write.table()
head(sleuth_significant, 20)
# The table shown above displays the top 20 significant genes 
# with a (Benjamini-Hochberg multiple testing corrected) q-value <= 0.05

# Resultados de isoformas
sleuth_results(so, 'conditionctl', test_type = "wt", which_model = "full", rename_cols = TRUE, show_all = TRUE)

# Resultados de genes
#sleuth_gene_table()

# Top 20 diferencialmente expressos
#head(table_lrt_res[, c(1, 4, 5, 16, 17)], n=20)


# Escrever uma tabela com resultados de sleuth
write.table(as.data.frame(head(sleuth_significant, 20)), sep = "\t", "sleuth_significant.txt", row.names = F)

##----------------------------------------------------------------
# Obten??o nome genes 
# get the gene names using biomaRt
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

write.table(as.data.frame(t2g), sep = ",", "sleuth_gene-names_WT.txt", row.names = F)


