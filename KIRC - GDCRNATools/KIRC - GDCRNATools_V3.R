##### KIRC - GDCRNATools #####
#####  Epitácio Farias  #####


####### Instalação dos Pacotes #######
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("GDCRNATools")

install.packages("gprofiler2")
install.packages("tidyverse")

####### Carregando os pacotes #######
library(GDCRNATools)
library(gprofiler2)
library(tidyverse)

####### Baixando os dados clínicos e de expressão (RNA-Seq e miRNA-Seq) do GDC #######
project <- 'TCGA-KIRC'
rnadir <- paste(project, 'RNAseq', sep='/')
mirdir <- paste(project, 'miRNAs', sep='/')

### Dado de RNAseq
gdcRNADownload(project.id     = 'TCGA-KIRC', 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = rnadir)

### Dado de miRNA maduro
gdcRNADownload(project.id     = 'TCGA-KIRC', 
               data.type      = 'miRNAs', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = mirdir)

### Dado Clínico
clinicaldir <- paste(project, 'Clinical', sep='/')
gdcClinicalDownload(project.id     = 'TCGA-KIRC', 
                    write.manifest = FALSE,
                    method         = 'gdc-client',
                    directory      = clinicaldir)


####### Organização dos dados #######

### Analisar metadado de RNAseq
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-KIRC',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)

### Filtrar amostras duplicadas em metadado de RNASeq 
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)

### Filtrar amostras de Tumor não primário e Tecido normal não sólido em metadado de RNAseq 
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)


### Analisar metadado de miRNAs
metaMatrix.MIR <- gdcParseMetadata(project.id = 'TCGA-KIRC',
                                   data.type  = 'miRNAs', 
                                   write.meta = FALSE)

### Filtrar amostras duplicadas em metadado de miRNAs
metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)

### Filtrar amostras de Tumor não primário e Tecido normal não sólido em metadado de miRNAs 
metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)


####### Mesclar dados brutos de contagem #######
### Mesclar dados de RNAseq 

## Baixando os dados do Xenabrowser
url <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.htseq_counts.tsv.gz"
destfile <- "kirc_counts.tsv.gz"
download.file(url, destfile)

## Extraindo os dados 
rnaCounts <- read_tsv(gzfile("kirc_counts.tsv.gz"))

## Alterando os nomes das observações, retirando a versão do Ensembl_ID 
rnaCounts$Ensembl_ID <- sub("\\..*", "", rnaCounts$Ensembl_ID)
rnaCounts <- column_to_rownames(rnaCounts, "Ensembl_ID")
rnaCounts <- as.data.frame(rnaCounts)

## Reverter a normalização log2(count + 1)
rnaCounts <- 2^(rnaCounts)-1

## Alterando os nomes das colunas, retirando o último dígito 
colnames(rnaCounts) <- substr(colnames(rnaCounts), 1,15)

## Buscando as colunas em comum entre os dados da metamatriz e dos dados de contagem
cols <- intersect(metaMatrix.RNA$sample, colnames(rnaCounts))
rnaCounts <- rnaCounts %>% dplyr::select(cols)

## Igualando o tamanho da metamatriz com o dado de contagem
metaMatrix.RNA <- metaMatrix.RNA[metaMatrix.RNA$sample %in% cols, ]

### Mesclar dados de miRNAs 
mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,
                         path      = mirdir, # the folder in which the data stored
                         organized = FALSE, # if the data are in separate folders
                         data.type = 'miRNAs')


### Normalização TMM e Transformação Voom
## Normalização dos dados de RNAseq 
rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)

## Normalização dos dados de miRNAs
mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)


####### Análise de Expressão Diferencial #######

### Converter para inteiro para rodar o DESeq2
rnaCounts[,1:ncol(rnaCounts)]=lapply(1:ncol(rnaCounts),function(x) {
  tryCatch({
    as.integer(rnaCounts[[x]])
  },warning = function(w) {
    rnaCounts[[x]]}
  )} )

### Utilizando o DESeq2
DEGAll_RNAs <- gdcDEAnalysis(counts     = rnaCounts, 
                             group      = metaMatrix.RNA$sample_type, 
                             comparison = 'PrimaryTumor-SolidTissueNormal', 
                             method     = 'DESeq2')


### Utilizando o DESeq2 para os dados de miRNAs
DEGAll_MIR <- gdcDEAnalysis(counts     = mirCounts, 
                            group      = metaMatrix.MIR$sample_type, 
                            comparison = 'PrimaryTumor-SolidTissueNormal', 
                            method     = 'DESeq2')


### Expressão diferencial de todos os dados
deALL <- gdcDEReport(deg = DEGAll_RNAs, gene.type = 'all')

### Expressão diferencial dos long non-coding
deLNC <- gdcDEReport(deg = DEGAll_RNAs, gene.type = 'long_non_coding')

### Expressão diferencial dos codificantes de proteínas
dePC <- gdcDEReport(deg = DEGAll_RNAs, gene.type = 'protein_coding')

### Expressão diferencial de todos os miRNAs
deALL_MIR <- gdcDEReport(deg = DEGAll_MIR, gene.type = 'all')

### Visualização das Expressões Diferenciais
## Volcanoplot para todos os RNAs

# Utilizando todos os RNAs diferencialmente expressos
cutoff <- sort(DEGAll$PValue)[10]
shrink.deseq.cut <- DEGAll %>% 
  mutate(TopGeneLabel=ifelse(PValue<=cutoff, symbol, ""))

# Plotando o volcano plot
ggplot(shrink.deseq.cut, aes(x = logFC, y= -log10(FDR))) + 
  geom_point(aes(colour=FDR < 0.01), pch=20, size=2) +
  labs(x="log Fold Change", y="-log10(FDR)") + 
  geom_label_repel(aes(label=TopGeneLabel), 
                   seed = 123,
                   max.time = 3,
                   max.iter = Inf,
                   size = 3,
                   box.padding = 2, 
                   max.overlaps = Inf)

## Volcanoplot para lncRNAs

# Utilizando todos os lncRNAs diferencialmente expressos
cutoff <- sort(deLNC$PValue)[10]
shrink.deseq_lncRNA.cut <- deLNC %>% 
  mutate(TopGeneLabel=ifelse(PValue<=cutoff, symbol, ""))

# Plotando o volcano plot
ggplot(shrink.deseq_lncRNA.cut, aes(x = logFC, y= -log10(FDR))) + 
  geom_point(aes(colour=FDR < 0.01), pch=20, size=2) +
  labs(x="log Fold Change", y="-log10(FDR)") + 
  geom_label_repel(aes(label=TopGeneLabel), 
                   seed = 123,
                   max.time = 3,
                   max.iter = Inf,
                   size = 3,
                   box.padding = 2, 
                   max.overlaps = Inf)

## Volcanoplot para miRNAs

# Utilizando os miRNAs diferencialmente expressos
cutoff <- sort(deALL_MIR$PValue)[10]
shrink.deseq_MIR.cut <- deALL_MIR %>% 
  mutate(TopGeneLabel=ifelse(PValue<=cutoff, rownames(deALL_MIR), ""))

# Plotando o volcano plot
ggplot(shrink.deseq_MIR.cut, aes(x = logFC, y= -log10(FDR))) + 
  geom_point(aes(colour=FDR < 0.01), pch=20, size=2) +
  labs(x="log Fold Change", y="-log10(FDR)") + 
  geom_label_repel(aes(label=TopGeneLabel), 
                   seed = 123,
                   max.time = 3,
                   max.iter = Inf,
                   size = 3,
                   box.padding = 2, 
                   max.overlaps = Inf)

## Barplot 
# RNAs
gdcBarPlot(deALL, angle = 45, data.type = 'RNAseq')

# lncRNAs
gdcBarPlot(deLNC, angle = 45, data.type = 'RNAseq')

# miRNAs
gdcBarPlot(deALL_MIR, angle = 45, data.type = 'miRNAS')

## Heatmap para RNAS
degName = rownames(deLNC)
gdcHeatmap(deg.id = degName, metadata = metaMatrix.RNA, rna.expr = rnaExpr)

## Heatmap para miRNAS
degName_MIR = rownames(deALL_MIR)
gdcHeatmap(deg.id = degName_MIR, metadata = metaMatrix.MIR, rna.expr = mirExpr)

####### Análise de Enriquecimento Funcional #######
## Acessando os nomes dos genes diferencialmente expressos para o enriquecimento funcional 
gene_information <- (deALL$symbol)

##Análise de enriquecimento funcional, baseada no organismo humano e com uma corredão baseada no FDR
gostres <- gost(gene_information, organism = "hsapiens",correction_method = "fdr")

## Plot do Enriquecimento funcional
#Iterativo
gostplot(gostres, capped = TRUE, interactive = TRUE)

#Estático
p <- gostplot(gostres, capped = FALSE, interactive = FALSE)

pp <- publish_gostplot(p, highlight_terms = c("MIRNA:hsa-miR-335-5p","WP:WP4018") , 
                       width = NA, height = NA, filename = NULL )
pp

### Visualizar mapas de vias em uma pagina online local
## Carregando pacote
library(pathview)

## Carregando informações para o shiny
deg <- deALL$logFC
names(deg) <- rownames(deALL)

pathways <- as.character(enrichOutput$Terms[enrichOutput$Category=='KEGG'])
pathways

shinyPathview(deg, pathways = pathways, directory = 'pathview')

####### Análise da rede de RNAs endógenos concorrentes (ceRNAs) #######
## Análise das redes de ceRNAs usando base de dados internas 
ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = 'spongeScan', 
                          pc.targets  = 'spongeScan', 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)


### Visualizando a rede de ceRNAs com o CYTOSCAPE
ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 & 
                        ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]

edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')

write.table(edges, file='edges.txt', sep='\t', quote=F)
write.table(nodes, file='nodes.txt', sep='\t', quote=F)

### Plot de Correlação
## Local
gdcCorPlot(gene1    = 'ENSG00000251165', 
           gene2    = 'ENSG00000091831', 
           rna.expr = rnaExpr, 
           metadata = metaMatrix.RNA)

## Em uma pagina online local
shinyCorPlot(gene1    = rownames(deLNC), 
             gene2    = rownames(dePC), 
             rna.expr = rnaExpr, 
             metadata = metaMatrix.RNA)


####### Análises Univariadas #######
### Análises dos riscos proporcionais de Cox
survOutput <- gdcSurvivalAnalysis(gene     = rownames(deALL), 
                                  method   = 'coxph', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA)

### Análise Kaplan-Meier
survOutput <- gdcSurvivalAnalysis(gene     = rownames(deALL), 
                                  method   = 'KM', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA, 
                                  sep      = 'median')

## Gráfico Kaplan-Meier
gdcKMPlot(gene     = 'ENSG00000130600',
          rna.expr = rnaExpr,
          metadata = metaMatrix.RNA,
          sep      = 'median')


## Gráfico Kplan-Meier em página local
shinyKMPlot(gene = rownames(deLNC), rna.expr = rnaExpr, 
            metadata = metaMatrix.RNA)
