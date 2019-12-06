## Análise em Nível de Transcrito Kallisto
## Utilizando Triplicatas
# Data: 30/11/2019
library(tximport)
library(biomaRt)
library(DESeq2)
library(readr)
library(dplyr)
library(rhdf5)

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


#### Parte IV - DESeq2
## Design com formula simples: sem avisos R; menor chance de erro.
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

# Agora, o objeto dds.Txi pode ser usado como aquele dds nos
# passos subsequentes de DESeq2.
head(dds.txi$replicate)

## Pre-filtering
# Filtrar por counts insignificantes.
keep <- rowSums(counts(dds.txi)) >= 10

# Renomear dds.txi para dds:
dds <- dds.txi[keep,]

# Observar
head(dds$replicate)

# Relevel factor para control como referencia
#reference <- 'control'
head(dds$condition, 9)

### Análise de Expressão Diferencial (DE)
# Objeto dds por DESeq2
dds <- DESeq(dds)

# A função results gera tabelas de resultados.
res <- results(dds)

# Visualizar
res

# Note que podemos especificar o coeficiente ou contraste 
# que queremos construir como uma tabela de resultados, usando:
res <- results(dds, contrast = c('condition', 'zika', 'control'))

# Visualizar
res

## Log fold change shrinkage for visualization and ranking¶
# Contração log fod change para visualização e ranqueamento.
# Shrinkage of effect size (LFC estimates)
resultsNames(dds)

# Para contrair (shrink) LFC passar objeto dds para função lfcShrink:
resLFC <- lfcShrink(dds, coef = 'condition_zika_vs_control', type = 'apeglm')

# Observar
resLFC

## Reordenando Resultados com p-values e adjusted p-values
# Ordenar os resultados da tabela por menor p value:
resOrdered <- res[order(res$pvalue), ]

# Aqui é possível exportar os resultados para um arquivo CSV.
# O arquivo de texto simples é criado com os resultados ordenados.
# Esse arquivo pode ser utilizado para observar os genes.

write.csv(as.data.frame(resOrdered), file="zika_vs_controls_resultsOrdered.csv")

# Summary
summary(res)

# Quantos adjusted p-values são menores que 0.1?
sum(res$padj < 0.1, na.rm = TRUE)

# FDR cutoff, alpha.
res05 <- results(dds, alpha=0.05)
summary(res05)

sum(res05$padj < 0.05, na.rm=TRUE)


### Independent Hypothesis Weighting
## Ponderação de Hipóteses Independentes
# Filtragem de p value: ponderar (weight) hipóteses para otimizar o poder.
# Está disponível no Bioconductor sob nome IHW.
resIHW <- results(dds, filterFun = ihw)
summary(resIHW)

# Quais menores que 0.1
sum(resIHW$padj < 0.1, na.rm=TRUE)

metadata(resIHW)$ihwResult

##### Parte V - Exploração de Resultados
## MA-plot
# A função plotMA mostra os log2 fold change atribuível a uma dada variável
# sobre a média de contagens normalizadas para todas as amostras no DESeqDataSet.
plotMA(res , ylim = c(-2, 2))

# Pontos em vermelho: se o valor p ajustado (adjusted p value) for menor que 0.1.

# Visualizando para log2 fold changes que foram contraídos (LFC)
# associado com mudanças log2 fold changes advindas de baixas contagens de genes
# sem requerimento de thresholds de filtragem arbitrários.
plotMA(resLFC, ylim=c(-2,2))

# Após utilizar plotMA pode-se utilizar a função identify para detectar interativamente
# o número de linhas de genes individuais ao clicar no plot.
# Pode-se então resgatar os IDs dos genes salvando os índices resultantes.
idx <- identify(res$baseMean, res$log2FoldChange)


### Alternative Shrinkage Estimators
## Alternative Shrinkage Estimators

# because we are interested in treated vs untreated, we set 'coef=2'
resNorm <- lfcShrink(dds, coef=2, type="normal")


## Lembrar de objeto LFC e Shrinkage


# Especificar o coeficiente pela ordem em que aparece em results(dds)
# O coeficiente usado em lfcShrink anterior (resNorm) foi "condition zika vs control"
# Porém, é possível especificar o coeficiente pela ordem em que aparece quando se usa resultsnames(dds):
resultsNames(dds)

# Usaremos o coeficiente como 2, pois é o que indica condition_zika_Vs_control.
# Nosso intersse se dá no contraste entre ambos:
# because we are interested in treated vs untreated, we set 'coef=2'
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")
resLFC <- lfcShrink(dds, coef = 'condition_zika_vs_control', type = 'apeglm')

# Agora, observar os plots juntos
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")


## Plot counts
# É útil examinar a contagem de reads para um único gene entre os grupos (control e zika).
# Existe a função plotCounts que pode fazer isso, a qual normaliza as contagens por profundidade
# de sequenciamento (sequencing depth) e adiciona uma pseudocontagem de 1/2 para permitir a plotagem
# em escala de log.
# Pode-se selecionar o gene de interesse a ser plotado por rowname ou por índice numérico.
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

# Customização com ggplot2
# Neste caso
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


# Outros genes: adicionar on ID (nome):
plotCounts(dds, gene='ENSG00000135845', intgroup="condition")

# Outros genes
e <- plotCounts(dds, gene='ENSG00000135845', intgroup="condition", returnData = T)
library("ggplot2")
ggplot(e, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

# Outros genes
f <- plotCounts(dds, gene='ENSG00000168300', intgroup="condition", returnData = T)
library("ggplot2")
ggplot(f, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

g <- plotCounts(dds, gene='ENSG00000254708', intgroup="condition", returnData = T)
library("ggplot2")
ggplot(g, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

## Mais informações na coluna Results
mcols(res)$description



############################## Outras observações ##############################

### DESeq2 oferece duas transformações para contabilizar dados, ambas estabilizam variâncias.
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
# Reordenada
# dds <- DESeq(dds)
res <- results(dds)
res
res <- results(dds, contrast = c('condition', 'zika', 'control'))
res

# Criar csv (pode ser usado em fgsea)
write.csv(as.data.frame(res), file = 'zika_vs_controls_results_res_GSEA.csv')

## Observação: há outra forma de criar
# Criar csv para fgsea para ctrstZikaxCTL:

#ctrstZikaxCTL <- as.data.frame(results(dds, contrast = c('condition','zika','control')))
#write.csv(ctrstZikaxCTL, 'zika_control_tximport_GSEA_2019.csv')

## Obejto res Reordenado por p-values e adjusted p-vaules:
resOrdered <- res[order(res$pvalue), ]

# Criar csv (pode ser usado em fgsea)
write.csv(as.data.frame(resOrdered), file="zika_vs_controls_results_reOrdered_GSEA.csv")

### Lembrar dos objetos:
# res
# resOrdered
# resLFC
# resNorm
# resAsh


#### Volcanoplot
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

