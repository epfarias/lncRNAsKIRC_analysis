# title: "R  Preprocessing and Download of Count Data - Breast Cancer TCGA"

# Installing and Loading Libraries            

if(!require("tidyverse")){install.packages("tidyverse")}
if(!require("readr")){install.packages("readr")}
if(!require("dplyr")){install.packages("dplyr")}

install.packages("dplyr")

# Downloading TCGA-BRCA clinical data from Xenabrowser

# Survival data
# https://xenabrowser.net/datapages/?dataset=TCGA-kirc.survival.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

url <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.survival.tsv"
destfile <- "kirc_survival.tsv.gz"
download.file(url, destfile)
kirc.survi <- read_tsv(gzfile("kirc_survival.tsv.gz"))
kirc.survi <- kirc.survi %>% 
  mutate(sample = str_replace_all(sample, "-", ".")) %>% 
  column_to_rownames("sample") %>% 
  rename(status = OS, obs.time = OS.time, patient_id = '_PATIENT')

kirc.survi <- as.data.frame(kirc.survi)

# Clinical data
# https://xenabrowser.net/datapages/?dataset=TCGA-BRCA.GDC_phenotype.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
url <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.GDC_phenotype.tsv.gz"
destfile <- "kirc_clinical.tsv.gz"
download.file(url, destfile)
kirc.clini <- read_tsv(gzfile("kirc_clinical.tsv.gz"))
kirc.clini <- kirc.clini %>%
  dplyr::select(c("submitter_id.samples","prior_malignancy.diagnoses","age_at_initial_pathologic_diagnosis", "gender.demographic",
                  "sample_type_id.samples", "pathologic_M", "pathologic_N", "pathologic_T", "ethnicity.demographic", "race.demographic")) %>% 
  dplyr::rename(sample = 'submitter_id.samples', 
                prior.dx = 'prior_malignancy.diagnoses', 
                age = 'age_at_initial_pathologic_diagnosis', 
                gender = 'gender.demographic',
                sample.type = 'sample_type_id.samples',
                metastasis = 'pathologic_M',
                neoplasm = 'pathologic_N',
                ajcc.stage = 'pathologic_T',
                ethnicity = 'ethnicity.demographic',
                race = 'race.demographic') %>% 
  dplyr::mutate(sample.type = str_replace_all(sample.type, "01", "TP") ) %>% 
  dplyr::mutate(sample.type = str_replace_all(sample.type, "11", "NT") ) %>% 
  dplyr::filter(sample.type %in% c("TP", "NT")) %>%  
  dplyr::mutate(sample = str_replace_all(sample, "-", ".")) %>% 
  dplyr::filter(sample %in% row.names(kirc.survi)) %>% 
  column_to_rownames(var = "sample") 

kirc.clini <- cbind(kirc.clini, kirc.survi[rownames(kirc.clini),])
kirc.clini$codes <- rownames(kirc.clini)


# Downloading TCGA-BRCA counts from Xenabrowser

# https://xenabrowser.net/datapages/?dataset=TCGA-BRCA.htseq_counts.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
url <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.htseq_counts.tsv.gz"
destfile <- "kirc_counts.tsv.gz"
download.file(url, destfile)
kirc.count <- read_tsv(gzfile(destfile))
kirc.count <- as.data.frame(kirc.count)
colnames(kirc.count) <- gsub("-", "\\.", colnames(kirc.count))
row.names(kirc.count) <- sub("\\..*", "", kirc.count$Ensembl_ID)
kirc.count$Ensembl_ID <- NULL
kirc.count <- 2^(kirc.count)-1
kirc.count <- round(kirc.count, digits = 0)

# Select anotation dataset

# If you haven't already installed devtools...
# Use devtools to install the package

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

install.packages("devtools")
devtools::install_github("stephenturner/annotables")
library(annotables)

kirc.annot <- grch38 %>%
  dplyr::filter(grch38$ensgene %in%  row.names(kirc.count)) %>%
  dplyr::select(ensgene, symbol, description)

# Filtering Counts and Clinical data

kirc.annot <- kirc.annot[!duplicated(kirc.annot$symbol), ]
kirc.count <- kirc.count[kirc.annot$ensgene,]
rownames(kirc.count) <- kirc.annot$symbol
kirc.clini <- kirc.clini[rownames(kirc.clini) %in% colnames(kirc.count), ]
kirc.clini2 <- kirc.clini2[rownames(kirc.clini2) %in% colnames(kirc.count), ]
kirc.count <- kirc.count[,colnames(kirc.count) %in% rownames(kirc.clini2)]

save(kirc.count, kirc.clini, kirc.clini2, kirc.annot, file="Data/brca_count.RData", compress=T)