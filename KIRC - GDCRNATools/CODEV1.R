if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GDCRNATools")
library(GDCRNATools)

project <- 'TCGA-KIRC'
rnadir <- paste(project, 'RNAseq', sep='/')
mirdir <- paste(project, 'miRNAs', sep='/')

### Download RNAseq data
gdcRNADownload(project.id     = 'TCGA-KIRC', 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               method = 'gdc-client', ## use gdc-client tool to download data
               directory      = rnadir)

### Download miRNAs data
gdcRNADownload(project.id     = 'TCGA-KIRC', 
               data.type      = 'miRNAs', 
               write.manifest = FALSE,
               method = 'gdc-client', ## use gdc-client tool to download data
               directory      = mirdir)