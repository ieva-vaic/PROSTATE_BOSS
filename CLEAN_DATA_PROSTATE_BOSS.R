#Prostate cancer project of gene expression and sPD(L) - clean data

#libraries + set directory
library(readxl)
library(tidyverse)
library(reshape2) 
library(rstatix) 
library(ggprism)
library(grid)
library(gtable)
library(glmnet)
library(pROC)
library(formattable)
library(gtsummary)
library(gt)
setwd("~/rprojects/OTHER DATA/BOSO data")

#load data################################################
#gene expression and clinical data
boss_data_ieva <- readxl::read_xlsx("Genu_raiskos_duomenys_is2018_v240729_nepersislinkes.xlsx")
#imunological and clinical data
boss_dataL_margarita <- readxl::read_xlsx("magraritos_data_nutrumpintas.xlsx", na = ".")  #"." to fill NA 

#choose matching samples
boss_data <- boss_data_ieva %>%
  filter(MARGARITOS_IMTIS_v2 == "YES")
dim(boss_data) #72 cases

#choose matching clinical data
DN_kodas <- boss_data$DN
boss_CLINICAL <- boss_dataL_margarita[c(boss_dataL_margarita$DN %in% DN_kodas), ]
dim(boss_CLINICAL) #72 cases

#merge datasets
boso_sujungtas <- left_join(
  boss_data, boss_CLINICAL, by = "DN"
)
dim(boso_sujungtas) #72 rows

#add row names
rownames(boso_sujungtas) <- boso_sujungtas$DN

#outliers########################################################################
#look at the dataset
PSA_hprt_values <- c("PCA3_PSA_HPRT"                                                                                                    
                     ,"PSMA_PSA_HPRT"                                                                                                    
                     ,"TMERG_PSA_HPRT"                                                                                                   
                     , "AR_FL_PSA_HPRT"                                                                                                   
                     , "HOX6_PSA_HPRT"   )

boso_sujungtas <- boso_sujungtas %>%
  mutate(clin_sign_m = case_when(
    `Gleason "3+3"=6; "3+4"=7; "4+3"=7.5; "4+4"=8; "5+3"=8.5` == 6 ~ "Clinically insignificant",
    `Gleason "3+3"=6; "3+4"=7; "4+3"=7.5; "4+4"=8; "5+3"=8.5` == 7 ~ "Clinically insignificant",
    `Gleason "3+3"=6; "3+4"=7; "4+3"=7.5; "4+4"=8; "5+3"=8.5` %in% c(7.5, 8) ~ "Clinically significant"
  ))
PSA_HPRT <- melt(boso_sujungtas, id.vars='clin_sign_m',  measure.vars=PSA_hprt_values) #be outliers
ggplot(PSA_HPRT) + 
  geom_boxplot(aes(x=clin_sign_m, y=value, color=variable)) +
  ylab(label = "Gene expression, normalized to HRPT1, divided by PSA expression") +
  xlab(label = "PCa clinical significance")  

#find outliers
#0) separate data
psa_hprt_out <- boso_sujungtas$`PSMA_PSA_HPRT`
#1)  Calculate quartiles 
Q1 <- quantile(psa_hprt_out, 0.25, na.rm = TRUE)
Q3 <- quantile(psa_hprt_out, 0.75, na.rm = TRUE)

# Calculate IQR
IQR <- Q3 - Q1

# Define threshold
threshold <- 3 #cia pasirenki kiek daug ismest nori

# Identify outliers
lower_bound <- Q1 - threshold * IQR
upper_bound <- Q3 + threshold * IQR

outliers2 <- psa_hprt_out[psa_hprt_out < lower_bound | psa_hprt_out > upper_bound]

# Find indices of outliers
outlier_indices2 <- which(psa_hprt_out < lower_bound | psa_hprt_out > upper_bound)
print(outlier_indices2) #3 randa
#find outlier names
print(boso_sujungtas[outlier_indices2,])

#clean the outliers: final dataset is 69 cases
boso_sujungtas_cleaned_PSA <- boso_sujungtas[-outlier_indices2,]
dim(boso_sujungtas_cleaned_PSA) #69 cases

#remove outlier (from clinical data / immunology) DN-61: final dataset is 68 cases
which(boso_sujungtas_cleaned_PSA$DN == "DN-061")
boso_sujungtas_cleaned_PSA <- boso_sujungtas_cleaned_PSA[-8,]
which(boso_sujungtas_cleaned_PSA$DN == "DN-061")
dim(boso_sujungtas_cleaned_PSA) #68 cases

#fix main clinical features########################################################
#stage##
boso_sujungtas_cleaned_PSA$`Pooperaine patologine stadija: T2=2;T3=3;T4=4` 
table(boso_sujungtas_cleaned_PSA$`Pooperacine patologine stadija`)

boso_sujungtas_cleaned_PSA <- boso_sujungtas_cleaned_PSA %>%
  mutate(stage_number = str_extract(`Pooperacine patologine stadija`, "[0-9]+"))

# Convert the new column to numeric if needed
boso_sujungtas_cleaned_PSA$stage_number <- as.numeric(boso_sujungtas_cleaned_PSA$stage_number)
table(boso_sujungtas_cleaned_PSA$stage_number) #50 vs 19

#grade##
table(boss_CLINICAL$`Gleason "3+3"=6; "3+4"=7; "4+3"=7.5; "4+4"=8; "5+3"=8.5`)

boso_sujungtas_cleaned_PSA <- boso_sujungtas_cleaned_PSA %>%
  mutate(grade_label = case_when(
    `Gleason "3+3"=6; "3+4"=7; "4+3"=7.5; "4+4"=8; "5+3"=8.5` == 6 ~ "grade 1",
    `Gleason "3+3"=6; "3+4"=7; "4+3"=7.5; "4+4"=8; "5+3"=8.5` == 7 ~ "grade 2",
    `Gleason "3+3"=6; "3+4"=7; "4+3"=7.5; "4+4"=8; "5+3"=8.5` %in% c(7.5, 8) ~ "grade 3"
  ))

table(boso_sujungtas_cleaned_PSA$grade_label, useNA = "a") #18, 42, 9

#clin signf##
boso_sujungtas_cleaned_PSA <- boso_sujungtas_cleaned_PSA %>%
  mutate(clinical_significance = case_when(
    `Gleason "3+3"=6; "3+4"=7; "4+3"=7.5; "4+4"=8; "5+3"=8.5` == 6 ~ "Clinically insignificant",
    `Gleason "3+3"=6; "3+4"=7; "4+3"=7.5; "4+4"=8; "5+3"=8.5` == 7 ~ "Clinically insignificant",
    `Gleason "3+3"=6; "3+4"=7; "4+3"=7.5; "4+4"=8; "5+3"=8.5` %in% c(7.5, 8) ~ "Clinically significant"
  ))
table(boso_sujungtas_cleaned_PSA$clinical_significance, useNA = "a") #60 vs 9

#save as RDS
saveRDS(boso_sujungtas_cleaned_PSA, "boso_sujungtas_cleaned_PSA.RDS")
#read rds fo chek####################################
boso_sujungtas_cleaned_PSA <- readRDS("boso_sujungtas_cleaned_PSA.RDS")
