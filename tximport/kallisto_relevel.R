## Análise em Nível de Transcrito Kallisto
## Utilizando Triplicatas
# Data: 01/12/2019

library(DESeq2)
library(readr)
library(dplyr)
library(rhdf5)
library(biomaRt)
library(readr)
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

# Relevel factor para control como referencia
#reference <- 'control'
samples$condition <- relevel(x = samples$condition,
                             ref = 'control')


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

#### DESeq2
## Design com formula simples:
dds.txi <- DESeqDataSetFromTximport(txi = txi.kallisto,
                                    colData = samples,
                                    design = ~replicate + condition)

# Agora, o objeto dds.Txi pode ser usado como aquele dds nos
# passos subsequentes de DESeq2.
dds.txi$replicate

# Filtrar por counts insignificantes.
keep <- rowSums(counts(dds.txi)) >= 10

# Renomear dds.txi para dds:
dds <- dds.txi[keep,]
rm(keep)

# Objeto dds por DESeq2
dds <- DESeq(dds)

# DESeq2 oferece duas transformações para contabilizar dados, ambas estabilizam variâncias.
# A. rlog: regularized-logarithm transformation (rlog); 
# B. VST: variance stabilizing transformation (vst);

## RLD
# Criar objeto rld, que transforma os dados para log2FC. 
# Scale/fator para minimizar diferenças entre amostras.
rld <- rlogTransformation(dds)

# rld retorna objeto DESeqTransform o qual é baseado na classe SummarizedExperiment.
head(assay(rld), 6)

# Os valores transformados não são contagens, sendo armazenados no slot assay.
# colData está ligado a dds e é acessível:
colData(rld)

## VST
# Criar objeto VST (recomendado para amostras n > 30); muito rápido.
vsd <- vst(dds)

# VST retorna objeto DESeqTransform o qual é baseado na classe SummarizedExperiment.
head(assay(vsd), 6)

# Os valores transformados não são contagens, sendo armazenados no slot assay.
# colData está ligado a dds e é acessível:
colData(vsd)


### Análise Exploratória de Dados
## Principal Component Analysis - PCA

# Gerar PCA.
# rld
pca <- plotPCA(rld, 
               ntop = nrow(counts(dds)),
               returnData = FALSE)


# Visualizar a PCA rld:
pca

# VST
pca.vst <- plotPCA(vsd,
                   ntop = nrow(counts(dds)),
                   returnData = FALSE)


# Visualizar a PCA vst:
pca.vst

# Extrair os resultados da análise de DE.
ctrstZikaxCTL <- as.data.frame(results(dds, 
                                       contrast = c('condition','zika','control')))

# Condition é a coluna/variável condition.



###----------------------- GSEA tximport - Arquivo 1 ---------------------

# Criar csv para fgsea para ctrstZikaxCTL:
#write.csv(ctrstZikaxCTL, 'zika_control_tximport_GSEA_2019.csv')

# Input para análise GSEA com fgsea (R, Bioconductor)

###----------------------------------------------------------------------


# Volcanoplot
with(as.data.frame(ctrstZikaxCTL[!(-log10(ctrstZikaxCTL$padj) == 0), ]), 
     plot(log2FoldChange,-log10(padj), 
          pch=16, 
          axes=T, 
          xlim = c(-6,6), 
          ylim = c(0,4),
          xlab = NA, 
          ylab = "-Log10(Pvalue-Adjusted)", main = "Febre zika vs Não Infectados \n(~replicate + condition)"
          
     )
)

with(subset(subset(as.data.frame(ctrstZikaxCTL), padj <= 0.05),
            log2FoldChange <= -1), 
     points(log2FoldChange,-log10(padj),
            pch=21, 
            col="black",bg = "#69B1B7"))

with(subset(subset(as.data.frame(ctrstZikaxCTL), padj <= 0.05),
            log2FoldChange >= 1), 
     points(log2FoldChange,-log10(padj),
            pch=21,
            col="black",bg = "tomato3"))

abline(h=1.3,col="green", lty = 2, cex= 3)
abline(v=1,col="green", lty = 2, cex= 3)
abline(v=-1,col="green", lty = 2, cex= 3)

