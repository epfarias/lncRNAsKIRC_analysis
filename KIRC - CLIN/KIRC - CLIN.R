### TCGA Clinical Data Analysis ####
# Epitácio Farias - 04/2022

# Project: TCGA-KIRC (Kidney Renal Clear Cell Carcinoma)
# Disease Type: Adenomas and Adenocarcinomas
# Cases: 537


## Install packages
pkg.bioconductor <- c("TCGAbiolinks", "TCGAWorkflow")
pkg.cran <- c("tidyverse", "skimr", "tableone", "survminer", "survival","finalfit", "DT", "data.table","rio")
# optional pkgs: "finalfit", "DT", "data.table"

#check if each package is on the local machine
#if it is not installed, it will be installed and loaded
pkg.check <- lapply(pkg.bioconductor, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, dependencies = TRUE)
    library(x, character.only = TRUE)}
})
pkg.check <- lapply(pkg.cran, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)}
})

rm(pkg.cran, pkg.bioconductor, pkg.check)
# devtools::install_github("BioinformaticsFMRP/TCGAWorkflow")


## 1. Data importing and visualizing ---------------------------
# https://www.bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html
# XML files includes all clinical, biospecimen, follow-ups and admin information
# indexed clinical data (XML subset) includes diagnoses, treatments, demographic and exposures information
# biospecimen data (get sample types code for normal tissue)

# uvm.biosp <- GDCquery_clinic(project = "TCGA-UVM", type = "Biospecimen", save.csv = FALSE)
# Fetch clinical data directly from the clinical XML files, using clinical.info parameter
# uvm.clinic <- GDCquery(
#   project = "TCGA-UVM",
#   file.type = "xml",
#   data.category = "Clinical",
# ) 
# GDCdownload(uvm.clinic)

kirc.clinic <- GDCquery_clinic(project = "TCGA-KIRC", type = "Clinical", save.csv = FALSE)

# uvm.patient <- GDCprepare_clinic(uvm.clinic, clinical.info = "patient")
# datatable(clinic.patient, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

# uvm.drug <- GDCprepare_clinic(uvm.clinic, clinical.info = "drug")
# datatable(clinic.drug, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

# uvm.radiation <- GDCprepare_clinic(uvm.clinic, clinical.info = "radiation")
# datatable(clinic.radiation, options = list(scrollX = TRUE,  keys = TRUE), rownames = FALSE)

# uvm.admin <- GDCprepare_clinic(uvm.clinic, clinical.info = "admin")
# datatable(clinic.admin, options = list(scrollX = TRUE, keys = TRUE), rownames = FALSE)

glimpse(kirc.clinic)
view(kirc.clinic)


## 2. Cleaning data ---------------------------

# Cleaning names 
# uvm.clinic<- janitor::clean_names(uvm.clinic)
names(kirc.clinic)

kirc.names <- names(kirc.clinic) %>%
  str_remove("treatments_")
# str_replace_all("_", ".")
# rm_accent()
names(kirc.clinic) <- kirc.names


# Cleaning variables 
kirc.clinic %>% missing_glimpse()
skim(kirc.clinic) 

logical.vars <- names(kirc.clinic %>% select_if(is.logical))
kirc.clinic <- kirc.clinic %>% 
  select(-all_of(logical.vars))

# Select variables based on NA count (> 50% complete is a good choice!).
# 
# NA_fifty <- dim(uvm.clinic)[1]/2
# NA_sum <- colSums(is.na(uvm.clinic))
# NA_sum <- as.data.frame(NA_sum)
# NA_sum <- tibble::rownames_to_column(NA_sum, "variables")
# NA_sum <- NA_sum %>%
#   filter(NA_sum < NA_fifty)
# 
# uvm_clinic <- uvm.clinic %>%
#   select(one_of(NA_sum$variables))

## Remove duplicate observations (patient_id or other id variable)
kirc.clinic <- kirc.clinic %>%
  distinct_at('submitter_id', .keep_all = TRUE)

## Remove numeric variables with unique observations (or <50% complete)
# need to maintain 'days_to_death' for survival anaysis?
kirc.clinic %>%
  select_if(is.numeric) %>%
  skim()

kirc.clinic <- kirc.clinic  %>%
  select(!c('days_to_diagnosis', 'days_to_birth', 'days_to_death', 'year_of_death'))

## Remove character variables with unique observations 
# unique(uvm.clinic$variable.name)
kirc.clinic %>%
  select_if(is.character) %>%
  skim()

kirc.clinic <- kirc.clinic  %>%
  select(!c('last_known_disease_status', 'prior_treatment', 'state', 
            'ajcc_staging_system_edition', 
            'classification_of_tumor', 'tumor_grade', 
            'progression_or_recurrence', 'alcohol_history', 
            'pharmaceutical_treatment_type', 'radiation_treatment_type'))

## Remove variables with similar information - check each one!
kirc.clinic$submitter_id == kirc.clinic$bcr_patient_barcode
kirc.clinic$age.at.index == as.integer(kirc.clinic$age_at_diagnosis/365)
kirc.clinic$tissue.or.organ.origin == kirc.clinic$site_of_resection_or_biopsy

kirc.clinic <- kirc.clinic  %>%
  select(!c('bcr_patient_barcode', 'age_at_diagnosis', 'site_of_resection_or_biopsy')) 

## Check and remove other variables
names.id <- names(kirc.clinic %>%
                    select(ends_with("_id")))
names.id <- names.id[-1]

kirc.clinic <- kirc.clinic %>% 
  select(-all_of(names.id))


## 3. Changing variables names ---------------------------
# Use a significative name, with snake_style or other. 
# names.clean <- names(uvm.clinic) %>% 
#    str_remove(".of") %>% 
#    str_remove(".or.therapy")
# names(uvm.clinic) <- names.clean

kirc.clinic <- kirc.clinic %>%
  rename(patient_id = 'submitter_id', age= "age_at_index")


## 4. Taming data ------------------------------------------
kirc.clinic <- kirc.clinic %>%
  mutate_if(is.character, as.factor) %>%
  mutate(patient_id = as.character(patient_id))

kirc.clinic$updated_datetime <- lubridate::as_date(kirc.clinic$updated_datetime)


## 5. Checking NA patterns -----------------------------
# Check distincts types of NAs: MCAR, MAR, MNAR
kirc.clinic  %>%
  missing_plot()


## 6. Checking numeric variables -----------------------------
# Check data distribution, unplausible values, outliers.
# Never delete an unusual value if the value is a possible one. 
# Deleting unusual values will bias your results and cause you to underestimate the variability in the observations.
# If the median and mean are similar, the distribution is likely roughly symmetrical. 
# Otherwise, it will be skewed to the right or to the left.

kirc.clinic %>% 
  select_if(is.numeric) %>%
  summary()

# uvm.clinic.num <- uvm.clinic %>% 
#    select_if(is.numeric)

# Histograms or density plots
ggplot(kirc.clinic, aes(age)) +
  geom_histogram(bins = 20, alpha = 0.5, color = "red")

ggplot(kirc.clinic, aes(year_of_birth)) +
  geom_histogram(bins = 20, alpha = 0.5, color = "red")

ggplot(kirc.clinic, aes(year_of_diagnosis)) +
  geom_bar(stat= "count", alpha = 0.5, color = "blue")

ggplot(kirc.clinic, aes(days_to_last_follow_up)) +
  geom_density( alpha = 0.5, color = "purple")

# hist <- lapply(names(uvm.clinic.num), function(i) {hist(uvm.clinic.num[,i],100,
#                                                         col="light blue",
#                                                         main=paste0("Histogran of ",i),
#                                                         xlab=i)})

# Boxplots 
# Inter quartil range (IQR) = Q3 — Q1
# Outliers = Q1 − 1.5 ∗ IQR < valor < Q3 + 1.5 ∗ IQR

ggplot(kirc.clinic, aes(x ='', y=age)) +
  geom_boxplot(width = .5) +
  geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(kirc.clinic$age)

ggplot(kirc.clinic, aes(x ='', y=days_to_last_follow_up)) +
  geom_boxplot(width = .5) +
  geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(kirc.clinic$days_to_last_follow_up)

# boxplot <- lapply(names(uvm.clinic.num), function(i) {boxplot(uvm.clinic.num[,i],
#                                                               main=paste0("Boxplot of ",i),
#                                                               xlab=i)})


## 7. Checking categorical variables --------------------------
# Check frequency, lables and levels 
# Cancer staging: https://www.cancer.gov/about-cancer/diagnosis-staging/staging
kirc.clinic %>%
  select_if(is.factor) %>%
  summary() 

# agregating levels

kirc.clinic <- kirc.clinic %>%
  mutate(ajcc_pathologic_t = fct_collapse(ajcc_pathologic_t, T1=c('T1', 'T1a', 'T1b'), T3=c('T3a', 'T3b')))


# recode levels
# uvm.clinic <- uvm.clinic %>%
#    mutate(ajcc.clinical.t = fct_recode(ajcc.clinical.t, T2 = T2a))

# drop levels
#kirc.clinic <- kirc.clinic %>%
  #mutate(primary_diagnosis = fct_recode(primary_diagnosis, NULL="Malignant melanoma, NOS"))

# uvm.clinic <- uvm.clinic %>%
#      mutate(gender = if_else(gender %in% c('male', 'female'), 1, 0))

kirc.clinic <- kirc.clinic  %>%
  select(!c(morphology, icd_10_code))

## Graphics

# categorical frequencies graphics
kirc.clinic %>% 
  count(ajcc_pathologic_stage) %>% 
  knitr::kable()

kirc.clinic %>% 
  ggplot(aes(x = tissue_or_organ_of_origin, fill = tissue_or_organ_of_origin)) + 
  geom_bar() +
  theme_bw(15) +
  labs(title = 'Age frequency by Tissue Site', x = "tissue site", y = "age") + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

kirc.clinic %>% 
  ggplot(aes(x = tissue_or_organ_of_origin, y = age, fill = tissue_or_organ_of_origin)) + 
  geom_boxplot() +
  theme_bw(15) +
  labs(title = 'Tissue Site', x = "tissue site", y = "age") + 
  facet_wrap(~ gender) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')



## 8. Saving dataset ------------------------------------
skim(kirc.clinic)
write_csv(kirc.clinic, "kirc.clinic.csv")


## 9. Summary table with tableone ------------------------------------
kirc.clinic2 <- kirc.clinic %>% 
  select(!where((is.character))) %>% 
  select(!(updated_datetime)) %>% 
  mutate(days_to_last_follow_up = as.integer(days_to_last_follow_up),
         year_of_diagnosis = as.integer(year_of_diagnosis),
         year_of_birth = as.integer(year_of_birth))

CreateTableOne(data=kirc.clinic2)

myVars <- names(kirc.clinic2)
catVars <- names(kirc.clinic2 %>% 
                   select_if(is.factor))
biomarkers <- c("days_to_last_follow_up", "year_of_diagnosis", "year_of_birth")

kirc_tab <- CreateTableOne(vars=myVars, data=kirc.clinic2, factorVars=catVars)
print(kirc_tab, nonnormal=biomarkers, showAllLevels=TRUE, formatOptions=list(big.mark = ","))

## 10. Survival analysis ------------------------------------
kirc.clinic3 <- kirc.clinic

# TCGABiolinks pkg (do not run)
# TCGAanalyze_survival(uvm.clinic, "gender")

# Survival pkg
# Dichotomize and change data labels
kirc.clinic3$vital_status <- as.numeric(ifelse(kirc.clinic3$vital_status == "Alive", 0, 1))
head(kirc.clinic3$vital_status)

# Fit survival data using the Kaplan-Meier method
# Time records survival time. Status indicates if patient is death (status=1) or was censored (status=0).
# A “+” after the time in the print out of km indicates censoring.
surv_object <- Surv(time = as.numeric(kirc.clinic3$days_to_last_follow_up), event = kirc.clinic3$vital_status)
surv_object 

fit <- survfit(surv_object ~ ajcc_pathologic_stage, data = kirc.clinic3)
fit
summary(fit)
# summary(fit1, times = c(1,30,60,90*(1:10)))
summary(fit)$table

ggsurvplot(fit, data = kirc.clinic3, pval = TRUE)
ggsurvplot(fit, data = kirc.clinic3,
           conf.int = TRUE,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           linetype = "strata",
           palette = c("#E7B800", "#2E9FDF"),
           legend = "bottom",
           legend.title = "Race",
           legend.labs = c("Not Reported", "White","Asian","Black oR African American"))

# Fit a Cox proportional hazards model # ERROR
fit.coxph <- coxph(surv_object ~ ajcc_pathologic_stage + prior_malignancy, data = kirc.clinic3)
ggforest(fit.coxph, data = kirc.clinic3)


## References ------------------------------------

# Packages documentation and vignettes.
# Flexible Imputation of Missing Data: https://stefvanbuuren.name/fimd/ 
# “Tame” data principles and the fivethirtyeight R package: https://rstudio-pubs-static.s3.amazonaws.com/396363_adaf67178eab4bd793bd9dd17dda70b3.html#today%E2%80%99s_focus
# Cancer staging: https://www.cancer.gov/about-cancer/diagnosis-staging/staging
