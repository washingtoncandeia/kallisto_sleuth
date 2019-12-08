## Análise em Nível de Transcrito Kallisto
## Utilizando Triplicatas
# Wald test p-value: condition zika vs control 
# Data: 08/12/2019
library(tximport)
library(apeglm)
library(biomaRt)
library(DESeq2)
library(readr)
library(dplyr)
library(rhdf5)
library(IHW)
library(GenomicAlignments)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(genefilter)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ReportingTools)
library(Glimma)
library(pcaExplorer)

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
pop <- c(rep(GBS, 3), 
         rep(CONTROL, 6), 
         rep(ZIKA, 3), 
         rep(GBS, 6), 
         rep(GBS_REC, 3), 
         rep(CHIKV_REC, 3), 
         rep(CONTROL, 3), 
         rep(GBS, 3), 
         rep(CONTROL, 3),
         rep(CHIKV, 3), 
         rep(ZIKA, 3), 
         rep(GBS, 3), 
         rep(GBS_REC, 9), 
         rep(CHIKV, 6), 
         rep(CONTROL, 6), 
         rep(CHIKV_REC, 3), 
         rep(ZIKA, 3), 
         rep(GBS_REC, 6), 
         rep(GBS, 3), 
         rep(CHIKV, 3), 
         rep(CHIKV_REC, 12), 
         rep(ZIKA, 3), 
         rep(GBS_REC, 9))


head(pop, 12)    # ZIKA, CHIKV, CHIKV_REC, GBS, GBS_REC, CONTROL
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


condition <- c(rep('gbs', 3), 
               rep('control', 6), 
               rep('zika', 3), 
               rep('gbs', 6), 
               rep('gbs_rec', 3), 
               rep('chikv_rec', 3), 
               rep('control', 3), 
               rep('gbs', 3), 
               rep('control', 3),
               rep('chikv', 3), 
               rep('zika', 3), 
               rep('gbs', 3), 
               rep('gbs_rec', 9), 
               rep('chikv', 6), 
               rep('control', 6), 
               rep('chikv_rec', 3), 
               rep('zika', 3), 
               rep('gbs_rec', 6), 
               rep('gbs', 3), 
               rep('chikv', 3), 
               rep('chikv_rec', 12), 
               rep('zika', 3), 
               rep('gbs_rec', 9))


## Parte II
# Aqui inicia-se a construção do data frame com todas as informações de amostras.


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


# Observações gerais:
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
write.table(samples_info, 'condition_zika_vs_control_GSEA.txt', sep = '\t')


# Criar um vetor nomeado apontando os arquivos de quantificação.
# Estes arquivos têm seus nomes anotados em uma tabela (samples.txt).
samples <- read.table('condition_zika_vs_control_GSEA.txt', header = TRUE, row.names = 1)
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

## Usando biomaRt para nomear transcritos e genes
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


### Parte III - Tximport Utlizando Arquivos kallisto
## Quantificação de Abundâncias de Transcritos com kallisto
## Análise de Expressão Diferencial (DE) com DESeq2.


# Estimativa de contagens a partir de kallisto,
# Usar ignoreTxVersion e ignoreAfterBar para que o data frame de IDs de transcritos
# e genes do Ensembl tenham ignorados as versões e barras |, respectivamente.
txi.kallisto <- tximport(files, 
                         type = 'kallisto',
                         tx2gene = t2g,
                         ignoreTxVersion = TRUE,
                         ignoreAfterBar = TRUE)


# Observações gerais
names(txi.kallisto)
head(txi.kallisto$abundance)

# Salvamento de um objeto R para uso posterior:
#dir.create(path = "./count_estimates/")
#save(txi.kallisto, file = "./count_estimates/txi_count_Trip_estimates.Rdata")

#### Parte IV - DESeq2
## Design com formula simples:
dds.txi <- DESeqDataSetFromTximport(txi = txi.kallisto,
                                    colData = samples,
                                    design = ~condition)

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
# Relevel como exemplo:
#dds$condition <- relevel(dds$condition, ref = "control")


### Análise de Expressão Diferencial (DE)
# Objeto dds por DESeq2
dds <- DESeq(dds)

# A função results gera tabelas de resultados.
res <- results(dds)

# Visualizar
res

# Note que podemos especificar o coeficiente ou contraste 
# que queremos construir como uma tabela de resultados, usando:
res <- results(dds, contrast = c('condition', 'control', 'zika'))

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

# Após utilizar plotMA pode-se utilizar a função identify para detectar interativamente
# o número de linhas de genes individuais ao clicar no plot.
# Pode-se então resgatar os IDs dos genes salvando os índices resultantes.
#idx <- identify(res$baseMean, res$log2FoldChange)

# Pontos em vermelho: se o valor p ajustado (adjusted p value) for menor que 0.1.

# Visualizando para log2 fold changes que foram contraídos (LFC)
# associado com mudanças log2 fold changes advindas de baixas contagens de genes
# sem requerimento de thresholds de filtragem arbitrários.
plotMA(resLFC, ylim=c(-2,2))

### Alternative Shrinkage Estimators
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

# Exporting only the results which pass an adjusted p value threshold 
# can be accomplished with the subset function, followed by the write.csv function.
resSig <- subset(resOrdered, padj < 0.1)


###### Parte VI
## Transformação e Visualização de Dados
# Transformações de contagem

# VST - Variance Stabilizing Transformation 
vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)

# Os valores transformados não são contagens, sendo armazenados no slot assay.
# colData está ligado a dds e é acessível:
colData(vsd)

# RLD - Regularized log Transformation 
rld <- rlog(dds, blind=FALSE)

# rld retorna objeto DESeqTransform o qual é baseado na classe SummarizedExperiment.
head(assay(rld), 6)

# Os valores transformados não são contagens, sendo armazenados no slot assay.
# colData está ligado a dds e é acessível:
colData(rld)

## Efeitos das Transformações na Variância
# Plot de desvio padrão dos dados transformados através das amostras,
# contra a média, usando shifting logarithm transformation
# fornece log2(n + 1)
ntd <- normTransform(dds)

## library(vsn)
# 1. Objeto ntd
meanSdPlot(assay(ntd))

# 2. Objeto vsd
meanSdPlot(assay(vsd))

# 3. Objeto rld
meanSdPlot(assay(rld))

## Qualidade de Dados por Clusterização e Visualização
# Heatmap da matriz de contagem 

#library(pheatmap)
select <- order(rowMeans(counts(dds,normalized = TRUE)),
                decreasing = TRUE)[1:50]

# Utilizando variável condition e replicates
df <- as.data.frame(colData(dds)[,c("condition","replicate")])  # Todas as linhas (genes) e variáveis (colunas condition e replicate)

# Pheatmap (condição e suas replicatas)
pheatmap(assay(ntd)[select,], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = FALSE, annotation_col = df)


# Utilizando variável condition e replicates
df2 <- as.data.frame(colData(dds)[,c("condition","pop")])  # Todas as linhas (genes) e variáveis (colunas condition e replicate)

# Pheatmap (condição e populações)
pheatmap(assay(ntd)[select,], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = FALSE, annotation_col = df2)


## VSD object
# vsd - condition, replicate (df)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


# vsd - condition, pop (df2)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df2)

## rld object
# rld - condition, replicate (df)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# rld - condition, pop (df2)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df2)


### Heatmap das Distâncias Amostra-Amostra
# O utro uso de dados transformados: sample clustering.
# Usando a função dist para transposição de matriz de contagem transformada.
sampleDists <- dist(t(assay(vsd)))

# library(RColorBrewer)
sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$run, sep=" - ")

colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


### PCA - Principal component plot das amostras
# O plot PCA está relacionado à matriz de distância e evidencia as amostras no plano de 2D
# abrangidas por seus primeiros componentes principais.
# Esse gráfico é útil para visualizar o efeito geral de covariantes experimentais (e batch effects).


## Usando VST
# vsd object
plotPCA(vsd, intgroup=c("condition", "run"))

# vsd object
plotPCA(vsd, intgroup=c("condition", "replicate"))


## Usando RLT
# rld object
plotPCA(rld, intgroup=c("condition", "run"))

# rld
plotPCA(rld, intgroup=c("condition", "replicate"))


## Utilizando ggplot2 com os dados
## ggplot2
# rsd object
pcaData <- plotPCA(vsd, intgroup=c("condition", "replicate"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()



# Extrair os resultados da análise de DE.
# Reordenada
# dds <- DESeq(dds)
res <- results(dds, contrast = c('condition', 'zika', 'control'))
res

# Criar csv (pode ser usado em fgsea)
write.csv(as.data.frame(res), file = 'zika_vs_controls_results_res_GSEA.csv')

## Observação: há outra forma de criar
# Criar csv para fgsea como variável ctrstZikaxCTL:
ctrstZikaxCTL <- as.data.frame(res)
# Agora, escrever arquivo .csv
#write.csv(ctrstZikaxCTL, file = 'zika_vs_controls_results_res_GSEA.csv')


## Obejto res Reordenado por p-values e adjusted p-vaules:
resOrdered <- res[order(res$pvalue), ]

# Criar csv (pode ser usado em fgsea)
write.csv(resOrdered, file = "zika_vs_controls_results_reOrdered_GSEA.csv")

### Lembrar dos objetos:
# res
# resOrdered
# resLFC
# resNorm
# resAsh

#### Volcanoplot
with(as.data.frame(ctrstZikaxCTL[!(-log10(res$padj) == 0), ]), 
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

