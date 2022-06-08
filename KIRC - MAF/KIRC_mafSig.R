########## Mutational Signature ######### 
# TCGA-KIRC gene signature with Maftools


## Installing packages
packages_bioconductor <- c("TCGAbiolinks","maftools","TCGAmutations","BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg19","SummarizedExperiment")
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
#BiocManager::install("PoisonAlien/TCGAmutations")

rm(packages_cran, packages_bioconductor, package.check)

#setwd()

#Download Mutational Data
dlbc.mutect.maf <- GDCquery_Maf("DLBC", pipelines = "mutect2")


# Number of mutations on muse:  6094
#dim(DLBC.muse.maf)[1]
# Number of mutations on mutect2:  6406
#dim(DLBC.mutect.maf)[1]
# Number of mutations on varscan2:  6280
#dim(DLBC.varscan2.maf)[1]
# Number of mutations on somaticsniper:  4969
#dim(DLBC.somaticsniper.maf)[1]
# Number of Patients: 37
#length(unique(substr(DLBC.mutect.maf$Tumor_Sample_Barcode,1,12)))

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

dlbc.mutect.maf_clin <- read.maf(maf = dlbc.mutect.maf,  
                                 verbose = T, 
                                 isTCGA = T, 
                                 removeDuplicatedVariants = F)

#Trinucleotide Matrix
dlbc.mutect.tnm = trinucleotideMatrix(maf = dlbc.mutect.maf_clin, prefix = '', add = TRUE,
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
plotApobecDiff(tnm = dlbc.tnm, maf = dlbc.mutect.maf_clin, pVal = 0.5)

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

dlbc.sign = estimateSignatures(mat = dlbc.tnm, nTry = 10, pConstant = 0.1, plotBestFitRes = T, parallel = 2)

#- Legacy - Mutational Signatures (v2 - March 2015):
# https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt
# https://cancer.sanger.ac.uk/signatures_v2/Signature_patterns.png
# https://cancer.sanger.ac.uk/signatures_v2/matrix.png


#- Single Base Substitution (SBS) - Mutational Signatures (v3.1 - June 2020)
# https://cancer.sanger.ac.uk/cosmic/signatures/SBS/index.tt

# Analysis with 4 gene signatures
dlbc.sig = extractSignatures(mat = dlbc.tnm, n = 4, pConstant = 0.1,  parallel = 2)

#Compate against original 30 signatures 
dlbc.og30.cosm = compareSignatures(nmfRes = dlbc.sig, sig_db = "legacy")

#library('pheatmap')
pheatmap::pheatmap(mat = dlbc.og30.cosm$cosine_similarities, 
                   cluster_rows = FALSE, 
                   angle_col = "45",
                   cellwidth = 20, cellheight = 20,
                   width = 7, height=4,
                   main = "Cosine similarity against validated signatures - Legacy")

#Compate against updated version3 60 signatures 
dlbc.v4.cosm = compareSignatures(nmfRes = dlbc.sig, sig_db = "SBS")

#library('pheatmap')
pheatmap::pheatmap(mat = dlbc.v4.cosm$cosine_similarities, 
                   cluster_rows = FALSE, 
                   angle_col = "45",
                   cellwidth = 20, cellheight = 20,
                   width = 7, height=4,                   
                   main = "Cosine similarity against validated signatures - SBS")

maftools::plotSignatures(nmfRes = dlbc.sig, title_size = 0.9, sig_db = "legacy")

maftools::plotSignatures(nmfRes = dlbc.sig, title_size = 0.9, sig_db = "SBS")

#Signatures can further be assigned to samples and enrichment analysis can be
#performd using signatureEnrichment funtion, which identifies mutations enriched
#in every signature identified.


dlbc.se = signatureEnrichment(maf = dlbc.mutect.maf_clin, sig_res = dlbc.sig)

#Above results can be visualzied similar to clinical enrichments.

plotEnrichmentResults(enrich_res = dlbc.se, pVal = 0.05)

#-----------------------------------------
# Workflow for Mutations on Firehose with hg19 ----
#-----------------------------------------

#Download Mutational Data
# https://github.com/PoisonAlien/TCGAmutations/        
dlbc.firehose <- TCGAmutations::tcga_load(study = "DLBC", source = "Firehose")

#Requires BSgenome object
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)

dlbc.firehose_clin <- read.maf(maf = dlbc.firehose@data, 
                                   clinicalData=dlbc.firehose@clinical.data, 
                                   verbose = T, 
                                   isTCGA = T, 
                                   removeDuplicatedVariants = F)

#Trinucleotide Matrix
dlbc.fh.tnm = trinucleotideMatrix(maf = dlbc.firehose_clin, prefix = 'chr', add = TRUE,
                               ref_genome = "BSgenome.Hsapiens.UCSC.hg19")        

#APOBEC Differentiation by Trinucleotide Matrix
# ---APOBEC related mutations are enriched in  2.083 % of samples (APOBEC enrichment score > 2 ;  1  of  48  samples
#plotApobecDiff(tnm = dlbc.fh.tnm, maf = dlbc.firehose_clin, pVal = 0.5)

#Run main function with maximum 10 signatures. 
library('NMF')

dlbc.sign = estimateSignatures(mat = dlbc.fh.tnm, nTry = 10, pConstant = 0.1, plotBestFitRes = T, parallel = 2)

#- Legacy - Mutational Signatures (v2 - March 2015):
# https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt
# https://cancer.sanger.ac.uk/signatures_v2/Signature_patterns.png
# https://cancer.sanger.ac.uk/signatures_v2/matrix.png


#- Single Base Substitution (SBS) - Mutational Signatures (v3.1 - June 2020)
# https://cancer.sanger.ac.uk/cosmic/signatures/SBS/index.tt

# Analysis with 4 gene signatures
dlbc.sig = extractSignatures(mat = dlbc.fh.tnm, n = 4, pConstant = 0.1,  parallel = 2)

#Compate against original 30 signatures 
dlbc.og30.cosm = compareSignatures(nmfRes = dlbc.sig, sig_db = "legacy")

#library('pheatmap')
pheatmap::pheatmap(mat = dlbc.og30.cosm$cosine_similarities, 
                   cluster_rows = FALSE, 
                   angle_col = "45",
                   cellwidth = 20, cellheight = 20,
                   width = 7, height=4,
                   main = "Cosine similarity against validated signatures - Legacy")

#Compate against updated version3 60 signatures 
dlbc.v4.cosm = compareSignatures(nmfRes = dlbc.sig, sig_db = "SBS")

#library('pheatmap')
pheatmap::pheatmap(mat = dlbc.v4.cosm$cosine_similarities, 
                   cluster_rows = FALSE, 
                   angle_col = "45",
                   cellwidth = 20, cellheight = 20,
                   width = 7, height=4,                   
                   main = "Cosine similarity against validated signatures - SBS")

maftools::plotSignatures(nmfRes = dlbc.sig, title_size = 0.9, sig_db = "legacy")

maftools::plotSignatures(nmfRes = dlbc.sig, title_size = 0.9, sig_db = "SBS")

#Signatures can further be assigned to samples and enrichment analysis can be
#performd using signatureEnrichment funtion, which identifies mutations enriched
#in every signature identified.


dlbc.se = signatureEnrichment(maf = dlbc.mutect.maf_clin, sig_res = dlbc.sig)

#Above results can be visualzied similar to clinical enrichments.

plotEnrichmentResults(enrich_res = dlbc.se, pVal = 0.05)


# References 
     
#1. maftools : Summarize, Analyze and Visualize MAF Files. 
#http://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html#910_mutational_signatures
     
#2. Alexandrov, L.B., et al., Signatures of mutational processes in human cancer. Nature, 2013. 500(7463): p.415-21.
#https://www.nature.com/articles/nature12477
     
#3. Roberts SA, Lawrence MS, Klimczak LJ, et al.An APOBEC Cytidine Deaminase Mutagenesis Pattern is Widespread in Human Cancers.
#Nature genetics. 2013;45(9):970-976. doi:10.1038/ng.2702. https://pubmed.ncbi.nlm.nih.gov/23852170/

#4. Signatures of Mutational Processes in Human Cancer
#https://cancer.sanger.ac.uk/cosmic/signatures

#5. Kidney renal clear cell: Signatures 1,2,5,9,13,17.
#https://cancer.sanger.ac.uk/signatures_v2/matrix.png

#sessionInfo()