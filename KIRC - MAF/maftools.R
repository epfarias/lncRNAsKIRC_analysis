############# KIRC_maftools ##################
# MAF files from Kidney Renal Clear Cell Carcinoma

# maftools: Summarize, Analyze and Visualize Mutation Anotated Files (MAF) Files
# URL: https://www.bioconductor.org/packages/release/bioc/html/maftools.html
# version 2.4.05


## Installing packages
packages_bioconductor <- c("TCGAbiolinks","maftools","BSgenome.Hsapiens.UCSC.hg38","SummarizedExperiment")
packages_cran <- c("DT", "tidyverse", "stringr", "data.table", "pheatmap","NMF")

#use this function to check if each package is on the local machine
#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed from Bioconductor and loaded
package.check <- lapply(packages_bioconductor, FUN = function(x) {
        if (!require(x, character.only = TRUE)) {
                BiocManager::install(x, dependencies = TRUE)
                library(x, character.only = TRUE)
        }
})
package.check <- lapply(packages_cran, FUN = function(x) {
        if (!require(x, character.only = TRUE)) {
                install.packages(x, dependencies = TRUE)
                library(x, character.only = TRUE)
        }
})

rm(packages_cran, packages_bioconductor, package.check)

setwd("~/lncRNAsKIRC_analysis/KIRC - MAF")

## Reading KIRC Maf files ---------------------------

# download MAF aligned against hg38
# it saves data inside a GDCdata and project name directory 
query <- GDCquery(
        project = "TCGA-KIRC", 
        data.category = "Simple Nucleotide Variation", 
        access = "open", 
        legacy = FALSE, 
        workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query)

#Remover duplicatas dentro do query 
query_results <- as.data.frame(query[[1]][[1]])
query_results <- query_results[!duplicated(query_results$cases),]
query[[1]][[1]] <- query_results

# GDC
maf <- GDCprepare(query)
sort(colnames(maf))

# MAF object contains main maf file, summarized data and any associated sample annotations
kirc.maf <- read.maf(maf = maf, useAll = T) 

# checking
getSampleSummary(kirc.maf) #  samples
getGeneSummary(kirc.maf) #  genes (hugo)
getClinicalData(kirc.maf) #  samples, no clinical data  
getFields(kirc.maf) #  variables 

# writes an output file
write.mafSummary(maf = kirc.maf, basename = 'kirc.maf')


## Reading clinical indexed data ------------------

# same as in data portal
clinical <- GDCquery_clinic(project = "TCGA-KIRC", type = "clinical", save.csv = FALSE)
sort(colnames(clinical))

colnames(clinical)[1] <- "Tumor_Sample_Barcode"
# clinical$Overall_Survival_Status <- 1
# clinical$Overall_Survival_Status[which(clinical$vital_status == "Dead")] <- 0
clinical$time <- clinical$days_to_death
clinical$time[is.na(clinical$days_to_death)] <- clinical$days_to_last_follow_up[is.na(clinical$days_to_death)]

# create object for survival analysis 
kirc.mafclin <- read.maf(maf = maf, clinicalData = clinical, isTCGA = T)


## Vizualizing -----------------------------------

# displays variants in each sample and variant types summarized by Variant_Classification
plotmafSummary(maf = kirc.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

# oncoplot for top ten mutated genes (costumize oncoplots!)
oncoplot(maf = kirc.maf, top = 10)

# transition and transversions
mafKIRC.titv <- titv(maf = kirc.maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = mafKIRC.titv)

# lollipop plot for VHL, which is one of the most frequent mutated gene in LAML
lollipopPlot(maf = kirc.maf, gene = 'VHL', AACol = 'HGVSp_Short', showMutationRate = TRUE)
lollipopPlot(maf = kirc.mafclin, gene = 'VHL', AACol = 'HGVSp_Short', showDomainLabel = FALSE)

# rainfall plots highlights hyper-mutated genomic regions
# Kataegis: genomic segments containing six or more consecutive mutations with an average inter-mutation distance of less than or equal to 1,00 bp
rainfallPlot(maf = kirc.maf, detectChangePoints = TRUE, pointSize = 0.6)

# mutation load
kirc.mutload <- tcgaCompare(maf = kirc.maf, cohortName = 'Example-KIRC', logscale = TRUE, capture_size = 50)

# plots Variant Allele Frequencies
# clonal genes usually have mean allele frequency around ~50% assuming pure sample
# it looks for column t_vaf containing vaf information, if the field name is different from t_vaf, we can manually specify it using argument vafCol
plotVaf(maf = kirc.maf, vafCol = 'Variants')


## Analysis -------------------------------------

# exclusive/co-occurance event analysis on top 10 mutated genes (pair-wise Fisher’s Exact test)
somaticInteractions(maf = kirc.maf, top = 10, pvalue = c(0.05, 0.1))

# detecting cancer driver genes based on positional clustering
# most of the variants in cancer causing genes are enriched at few specific loci (aka hot-spots)
kirc.sig <- oncodrive(maf = kirc.maf, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = kirc.sig, fdrCutOff = 0.1, useFraction = TRUE) # AACol = 'Amino_acids'?

# pfam domain
kirc.pfam = pfamDomains(maf = kirc.maf, AACol = 'Protein_Change', top = 10) # idem


## Reading gistic or CNV -------------------------

all.lesions <- system.file("extdata", "all_lesions.conf_99.txt", package = "maftools")
amp.genes <- system.file("extdata", "amp_genes.conf_99.txt", package = "maftools")
del.genes <- system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
scores.gis <- system.file("extdata", "scores.gistic", package = "maftools")

kirc.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes,
                         gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = TRUE)

# GISTIC object
kirc.gistic

# checking 
getSampleSummary(kirc.gistic) # 191 samples (Tumor_Sample_Barcode)
getGeneSummary(kirc.gistic) # 2622 genes (hugo)
getCytobandSummary(kirc.gistic) # 16 cytobands

write.GisticSummary(gistic = kirc.gistic, basename = 'kirc.gistic')


## Vizualizing Gistic ---------------------------------

# genome plot
gisticChromPlot(gistic = kirc.gistic, markBands = "all")

# bubble plot
gisticBubblePlot(gistic = kirc.gistic)

# oncoplot 
gisticOncoPlot(gistic = kirc.gistic, clinicalData = getClinicalData(x = kirc.maf), 
               clinicalFeatures = 'FAB_classification', sortByAnnotation = TRUE, top = 10) # error 

gisticOncoPlot(gistic = kirc.gistic, top = 10)

## Survival Analysis ------------------------------

# it requires input data with Tumor_Sample_Barcode, binary event (1/0) and time to event.
mafSurvival(maf = kirc.mafclin, genes = 'VHL', 
            time = 'days_to_last_follow_up', Status = 'vital_status', isTCGA = TRUE)

# identify a set of genes (of size 2) to predict poor prognostic groups
prog_geneset <- survGroup(maf = kirc.mafclin, top = 10, geneSetSize = 2, 
                          time = "days_to_last_follow_up", Status = "vital_status", verbose = FALSE)

## Clinical enrichment analysis ------------------

aacs.enrich <- clinicalEnrichment(maf = kirc.mafclin, clinicalFeature = 'tumor_grade')

# identify mutations associated with clinicalFeature
# results are returned as a list. Significant associations p-value < 0.05
aacs.enrich$groupwise_comparision[p_value < 0.05] # it takes too long!!

plotEnrichmentResults(enrich_res = aacs.enrich, pVal = 0.05)


## Oncogenic Signaling Pathways --------------------

OncogenicPathways(maf = kirc.maf)

PlotOncogenicPathways(maf = kirc.maf, pathways = "RTK-RAS")
# tumor suppressor genes (red) and oncogenes (blue)


## Mutational signatures --------------------------
#Gene Signature of Kidney Cancer (TCGA-KIRC, Diffuse Large B-Cell Lymphoma) with Maftools

#- References: 
        
#1. maftools : Summarize, Analyze and Visualize MAF Files. 
#http://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html#910_mutational_signatures

#2. Alexandrov, L.B., et al., Signatures of mutational processes in human cancer.
#Nature, 2013. 500(7463): p. 415-21. 
#https://www.nature.com/articles/nature12477

#3. Roberts SA, Lawrence MS, Klimczak LJ, et al.
#An APOBEC Cytidine Deaminase Mutagenesis Pattern is Widespread in Human Cancers.
#Nature genetics. 2013;45(9):970-976. doi:10.1038/ng.2702. https://pubmed.ncbi.nlm.nih.gov/23852170/
        
#4. Signatures of Mutational Processes in Human Cancer
#https://cancer.sanger.ac.uk/cosmic/signatures

#5. Kidney renal clear cell: Signatures 1,2,5,9,13,17.
#https://cancer.sanger.ac.uk/signatures_v2/matrix.png

#Download Mutational Data

kirc.mutect.maf <- GDCquery_Maf("KIRC", pipelines = "mutect2")

# Number of mutations on muse:  3180
#dim(laml.muse.maf)[1]
# Number of mutations on mutect2:  9905
#dim(laml.mutect.maf)[1]
# Number of mutations on varscan2:  4118
#dim(laml.varscan2.maf)[1]
# Number of mutations on somaticsniper:  2749
#dim(laml.somaticsniper.maf)[1]
# Number of Patients: 143
length(unique(substr(kirc.mutect.maf$Tumor_Sample_Barcode,1,12)))

#We select the mutect2 pipeline, since it has the larger number of variants.

#?Every cancer, as it progresses leaves a signature characterized by specific
#pattern of nucleotide substitutions. Alexandrov et.al have shown such 
#mutational signatures, derived from over 7000 cancer samples 5. Such signatures
#can be extracted by decomposing matrix of nucleotide substitutions, classified
#into 96 substitution classes based on immediate bases surrounding the mutated
#base. Extracted signatures can also be compared to those validated signatures.?

#?First step in signature analysis is to obtain the adjacent bases surrounding
#the mutated base and form a mutation matrix. NOTE: Earlier versions of maftools
#required a fasta file as an input. But starting from 1.8.0, BSgenome objects
#are used for faster sequence extraction.? [1]

#Requires BSgenome object
library(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)

kirc.mutect.maf_clin <- read.maf(maf = kirc.mutect.maf, 
                                 clinicalData=clinical, 
                                 verbose = T, 
                                 isTCGA = T, 
                                 removeDuplicatedVariants = F)

#Trinucleotide Matrix
kirc.tnm = trinucleotideMatrix(maf = kirc.mutect.maf_clin, prefix = '',
                               ref_genome = "BSgenome.Hsapiens.UCSC.hg38")

#In humans/mammals the APOBEC help protect from viral infections. The APOBEC
#enzymes, when misregulated, are a major source of mutation in numerous cancer
#types.

#?We can also analyze the differences in mutational patterns between APOBEC
#enriched and non-APOBEC enriched samples. plotApobecDiff is a function which
#takes APOBEC enrichment scores estimated by trinucleotideMatrix and classifies
#samples into APOBEC enriched and non-APOBEC enriched. Once stratified, it
#compares these two groups to identify differentially altered genes.?[1]

#?Note that, LAML with no APOBEC enrichments, is not an ideal cohort for this
#sort of analysis and hence below plot is only for demonstration purpose.?[1]

#APOBEC Differentiation by Trinucleotide Matrix
plotApobecDiff(tnm = kirc.tnm, maf = kirc.mutect.maf_clin, pVal = 0.2)

#Signature analysis includes following steps.

#1. estimateSignatures - which runs NMF on a range of values and measures the
#goodness of fit - in terms of Cophenetic correlation.

#2. plotCophenetic - which draws an elblow plot and helps you to decide
#optimal number of signatures. Best possible signature is the value at which
#Cophenetic correlation drops significantly.

#3. extractSignatures - uses non-negative matrix factorization to decompose
#the matrix into n signatures. n is chosen based on the above two steps. 
#In case if you already have a good estimate of n, you can skip above two steps.

#4. compareSignatures - extracted signatures from above step can be compared 
#to known signatures11 from COSMIC database, and cosine similarity is 
#calculated to identify best match.

#5. plotSignatures - plots signatures

par(mar = c(2, 2, 2, 1))
plot(NA, xlim = c(1, 10), ylim = c(0, 30), frame.plot = FALSE, axes = FALSE, xlab = NA, ylab = NA)
rect(xleft = 3, ybottom = 28, xright = 7, ytop = 30, col = grDevices::adjustcolor("gray70", alpha.f = 0.6), lwd = 1.2, border = "maroon")
text(x = 5, y = 29, labels = "MAF", font = 2)
arrows(x0 = 5, y0 = 28, x1 = 5, y1 = 26, length = 0.1, lwd = 2)
text(x = 5, y = 25, labels = "trinucleotideMatrix()", font = 3)
arrows(x0 = 5, y0 = 24, x1 = 5, y1 = 21, length = 0.1, lwd = 2)
text(x = 5, y = 20, labels = "estimateSignatures()", font = 3)
arrows(x0 = 5, y0 = 19, x1 = 5, y1 = 16, length = 0.1, lwd = 2)
text(x = 5, y = 15, labels = "plotCophenetic()", font = 3)
arrows(x0 = 5, y0 = 14, x1 = 5, y1 = 11, length = 0.1, lwd = 2)
text(x = 5, y = 10, labels = "extractSignatures()", font = 3)
arrows(x0 = 5, y0 = 9, x1 = 5, y1 = 6, length = 0.1, lwd = 2)
text(x = 5, y = 5, labels = "compareSignatures()", font = 3)
arrows(x0 = 5, y0 = 4, x1 = 5, y1 = 1, length = 0.1, lwd = 2)
text(x = 5, y = 0, labels = "plotSignatures()", font = 3)

#?Draw elbow plot to visualize and decide optimal number of signatures from
#above results.?

#?Best possible value is the one at which the correlation value on the y-axis
#drops significantly. In this case it appears to be at n = 3. LAML is not an
#ideal example for signature analysis with its low mutation rate, but for
#solid tumors with higher mutation burden one could expect more signatures,
#provided sufficient number of samples.?


#Run main function with maximum 10 signatures. 

library('NMF')
kirc.sign = estimateSignatures(mat = kirc.tnm, nTry = 2, pConstant = 0.1, plotBestFitRes = T, parallel = 2)

plotCophenetic(res = kirc.sign)
#- Legacy - Mutational Signatures (v2 - March 2015):
#        https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt
#https://cancer.sanger.ac.uk/signatures_v2/Signature_patterns.png
#https://cancer.sanger.ac.uk/signatures_v2/matrix.png


#- Single Base Substitution (SBS) - Mutational Signatures (v3.1 - June 2020)
#https://cancer.sanger.ac.uk/cosmic/signatures/SBS/index.tt

# Analysis with 4 gene signatures
kirc.sig = extractSignatures(mat = kirc.tnm, n = 4, pConstant = 0.1,  parallel = 2)

#Compate against original 30 signatures 
kirc.og30.cosm = compareSignatures(nmfRes = kirc.sig, sig_db = "legacy")

library('pheatmap')
pheatmap::pheatmap(mat = kirc.og30.cosm$cosine_similarities, 
                   cluster_rows = FALSE, 
                   angle_col = "45",
                   cellwidth = 20, cellheight = 20,
                   width = 7, height=4,
                   main = "Cosine similarity against validated signatures - Legacy")

#Compate against updated version3 60 signatures 
kirc.v4.cosm = compareSignatures(nmfRes = kirc.sig, sig_db = "SBS")

#library('pheatmap')
pheatmap::pheatmap(mat = kirc.v4.cosm$cosine_similarities, 
                   cluster_rows = FALSE, 
                   angle_col = "45",
                   cellwidth = 20, cellheight = 20,
                   width = 7, height=4,                   
                   main = "Cosine similarity against validated signatures - SBS")

maftools::plotSignatures(nmfRes = kirc.sig, title_size = 0.9, sig_db = "legacy")

maftools::plotSignatures(nmfRes = kirc.sig, title_size = 0.9, sig_db = "SBS")

#Signatures can further be assigned to samples and enrichment analysis can be
#performd using signatureEnrichment funtion, which identifies mutations enriched
#in every signature identified.


kirc.se = signatureEnrichment(maf = kirc.mutect.maf_clin, sig_res = kirc.sig)

#Above results can be visualzied similar to clinical enrichments.

plotEnrichmentResults(enrich_res = kirc.se, pVal = 0.05)
-----------------------------------------
-----------------------------------------
     
## References
     
# citation("maftools")
# Mayakonda A, Lin DC, Assenov Y, Plass C, Koeffler HP. 2018. Maftools: efficient and comprehensive analysis of somatic variants in cancer. Genome Resarch PMID: 30341162
     
# about MAF files: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
# about download data: https://www.bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/query.html
# about costumize oncoplots: http://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/oncoplots.html
# about validated signatures: https://cancer.sanger.ac.uk/cosmic/signatures
     
# sessionInfo()
