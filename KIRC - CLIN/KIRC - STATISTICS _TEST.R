######### TCGA Clinical Data Analysis #############
# checking variable relations to "vital_status"
# using Tidyverse whenever possible


install.packages("tidyverse", dependencies = TRUE)
install.packages("skimr", dependencies = TRUE)
install.packages("finalfit", dependencies = TRUE)
library(tidyverse)
library(skimr)
library(finalfit)

setwd()


## 1. Importing data ---------------------

kirc.clinic<- read_csv("kirc.clinic.csv")

## 2. Taming data -------------------------
# use lubridate for dates
kirc.clinic <- kirc.clinic %>%
  mutate_if(is.character, as.factor) %>%
  mutate(patient_id = as.character(patient_id),
         age = as.integer(age),
         year_of_diagnosis = as.integer(year_of_diagnosis))

kirc.clinic %>%
  select_if(is.factor) %>%
  skim()

kirc.clinic <- kirc.clinic  %>%
  select(!c("tissue_or_organ_of_origin","icd_10_code"))

glimpse(kirc.clinic)
View(kirc.clinic)

## 3. The dependent variable --------------------
# Check the number of levels. If greater than 2, run a ordinal logistic regression
# The independent variables are also called predicted or explanatory 
table(kirc.clinic$vital_status, exclude = NULL)
table(kirc.clinic$ajcc_pathologic_stage, exclude = NULL)

## 4. Categorical variables (chi-square test) --------------------
# Examine the relationship between chategorical variables and the dependent one, using cross tabulation.

# In the function prop.table(), using margin = 1 the table presents the row percentages, i.e. the frequencies of each column (vital_status) in each row (explanatory variable). If we entered margin = 2, this would display the inverse, i.e. the frequencies of each row in each column.

# ToDo: table with all 
t_status <- table(kirc.clinic$ajcc_pathologic_stage,kirc.clinic$vital_status, exclude = NA)
t_status <- addmargins(round(100*prop.table(t_status)))
t_status
chisq.test(x = kirc.clinic$ajcc_pathologic_stage, y = kirc.clinic$vital_status)


# summarise categorical variables by a categorical variable
explanatory_char <- kirc.clinic%>% 
  select_if(is.factor) %>%
  select(-vital_status) %>%
  names
dependent <- 'vital_status'

table_char <- kirc.clinic%>%
  summary_factorlist(dependent, explanatory_char, p=TRUE, 
                     na_include = FALSE, add_dependent_label=TRUE)
# To include NA into your statistic test: na_to_p = TRUE

table_char
# knitr::kable(table_char, row.names=FALSE, align=c("l", "l", "r", "r", "r"))

warnings()
# Droping levels with narrow distributions -> check warnings ()
# Group some levels or drop one (NULL = 'level') when grouping is not possible  


## 6. Numeric variables (t.test) ---------------
# Examine the relationship between your candidate predictor variables
# For continuous variables, use pairwise correlations and scatterplot matrices

# kirc.clinic_x <- kirc.clinic%>%
#      group_by(vital_status) %>%
#      select_if(is.numeric) %>%
#      summarise_at(vars(age:over_surv_mth), mean, na.rm = TRUE) %>%
#      group_by(vital_status)

# ggplot(kirc.clinic, aes(age, fill= vital_status)) +
#   geom_histogram(bins = 15, position = "dodge")
# t.test(kirc.clinic$age ~ kirc.clinic$vital_status)

# ggplot(kirc.clinic, aes(x=vital_status, y=disease_free_mth)) +
#   geom_boxplot(width = .5) +
#   geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
# t.test(kirc.clinic$disease_free_mth ~ kirc.clinic$vital_status)


# summarise numeric variables by a categorical variable
explanatory_num <- kirc.clinic%>% 
  select_if(is.numeric) %>%
  select(-year_of_diagnosis) %>%
  names
dependent <- 'vital_status'

table_fit_num <- kirc.clinic%>%
  summary_factorlist(dependent, explanatory_num, p=TRUE, 
                     na_include = TRUE, add_dependent_label=TRUE)
# to include NA into your statistic test: na_to_p = TRUE

table_fit_num
# knitr::kable(table_fit, row.names=FALSE, align=c("l", "l", "r", "r", "r"))

warnings()

# Correlation Matrix
# Use Pearson's (normal distribution) and Spearman (not-normal) correlations
# For continuous variables with no discontinuations and no obvious outliers
# Check collinearity (a strong linear relationship between predictors that causes problems in estimation of model parameters)
corr_num <- kirc.clinic%>%
  select_if(is.numeric) %>%
  drop_na()

# Check the correlation between variables to exclude the higly correlated
cor_matrix <- cor(corr_num, method = "spearman")
cor_matrix <- round(cor_matrix, 2)
cor_matrix

# To visualize
pairs(~ days_to_last_follow_up + year_of_diagnosis + age + year_of_birth, data=kirc.clinic)
# Higlhy collinears: over_surv_mt & disease_free_mth 
