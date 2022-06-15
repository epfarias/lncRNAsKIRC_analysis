##### KIRC - GDCRNATools #####
#####  Epitácio Farias  #####


####### Instalação do Pacote "GDCRNATools" #######
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("GDCRNATools")


####### Carregando o pacote #######
library(GDCRNATools)


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

library(tidyverse)

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
rnaCounts <- rnaCounts %>% select(cols)

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
### Utilizando o limma
DEGAll_limma <- gdcDEAnalysis(counts     = rnaCounts, 
                              group      = metaMatrix.RNA$sample_type, 
                              comparison = 'PrimaryTumor-SolidTissueNormal', 
                              method     = 'limma')

### Utilizando o DESeq2
DEGAll_DESeq <- gdcDEAnalysis(counts     = rnaCounts, 
                              group      = metaMatrix.RNA$sample_type, 
                              comparison = 'PrimaryTumor-SolidTissueNormal', 
                              method     = 'DESeq2')

### Carregando os dados de expressão diferencial
data(DEGAll)

### Expressão diferencial de todos os dados
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all')

### Expressão diferencial dos long non codign
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding')

### Expressão diferencial dos codificantes de proteínas
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding')

### Expressão diferencial dos miRNAS
demiR <- gdcDEReport(deg = DEGAll, gene.type = 'miRNAS')


### Visualização das Expressões Diferenciais
## Volcanoplot 
gdcVolcanoPlot(DEGAll)

## Barplot 
gdcBarPlot(deg = deALL, angle = 45, data.type = 'RNASeq')

## Heatmap
degName = rownames(deALL)
gdcHeatmap(deg.id = degName, metadata = metaMatrix.RNA, rna.expr = rnaExpr)

####### Análise de Enriquecimento Funcional #######
### Acessando e carregando os dados de enriquecimento
enrichOutput <- gdcEnrichAnalysis(gene = rownames(deALL), simplify = TRUE)
data(enrichOutput)

### Gráfico de barras
gdcEnrichPlot(enrichOutput, type = 'bar', category = 'GO', num.terms = 10)

### Gráfico de bolhas
gdcEnrichPlot(enrichOutput, type = 'bubble', category = 'GO', num.terms = 10)


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
                          lnc.targets = 'starBase', 
                          pc.targets  = 'starBase', 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)


### Análise das redes de ceRNAs usando base de dados providas pelo usuário
## carregando interações miRNA-lncRNA
data(lncTarget)

## carregando interações miRNA-mRNA
data(pcTarget)
pcTarget[1:3]

ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = lncTarget, 
                          pc.targets  = pcTarget, 
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
gdcKMPlot(gene     = 'ENSG00000136193',
          rna.expr = rnaExpr,
          metadata = metaMatrix.RNA,
          sep      = 'median')

## Gráfico Kplan-Meier em página local
shinyKMPlot(gene = rownames(deALL), rna.expr = rnaExpr, 
            metadata = metaMatrix.RNA)
