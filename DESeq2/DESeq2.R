##----------------------------------------
# Tese Cap. 2
# Febre zika x Controles
# Data: 01/12/2019
# Washington C. Araujo
# Análise a partir de tximport
##----------------------------------------
library(tximport)
library(DESeq2)
library(readr)
library(dplyr)
library(rhdf5)

# Iniciar a análise de tximport considerando dds.txi como o dds de DESeq2.
### Usando tximport para kallisto
load('./count_estimates/txi_count_Trip_estimates.Rdata')

# Observar o objeto txi.kallisto
head(txi.kallisto)
names(txi.kallisto)


#### DESeq2
## Design com formula simples:
dds.txi <- DESeqDataSetFromTximport(txi = txi.kallisto,
                                    colData = samples,
                                    design = ~condition)

## Um objeto em DESeq2 seria:
# 8. Gerar entradas para DESeq.
#dds <- DESeqDataSetFromMatrix(countData = countdata, 
#                              colData = coldata,
#                              design = ~condition)


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
write.csv(ctrstZikaxCTL, 'zika_control_tximport_GSEA_2019.csv')

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
          ylab = "-Log10(Pvalue-Adjusted)", main = "Febre zika vs Não Infectados"
                                                                            
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

