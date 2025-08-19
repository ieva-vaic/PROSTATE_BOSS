#Prostate cancer project of gene expression and sPD(L) - HEATMAP

#libraries + set directory
library(Hmisc) 
library(pROC)
library(survival)
library(survminer)
library(glmnet)
library(ggplotify)
library(formattable)
library(gtsummary)
library(gt)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)
setwd("~/rprojects/OTHER DATA/BOSO data")

#load data
boso_sujungtas_cleaned_PSA <- readRDS("boso_sujungtas_cleaned_PSA.RDS")
PSA_hprt_values <- c("PCA3_PSA_HPRT" ,"PSMA_PSA_HPRT" ,"TMERG_PSA_HPRT" ,
                     "AR_FL_PSA_HPRT" , "HOX6_PSA_HPRT"   )


PSA_hprt_values3 <- c("PCA3_PSA_HPRT"                                                                                                    
                      ,"PSMA_PSA_HPRT"                                                                                                    
                      , "AR_FL_PSA_HPRT"                                                                                                   
)

main_biomarkers <- c("PCA3_PSA_HPRT"                                                                                                    
                     ,"PSMA_PSA_HPRT"                                                                                                    
                     , "AR_FL_PSA_HPRT"
                     , "I sPDL1, pg/ml"
                     , "I sPD1, pg/ml"
)

boso_sujungtas_cleaned_PSA$Stage <- factor(boso_sujungtas_cleaned_PSA$Stadija_skaicius)
levels(boso_sujungtas_cleaned_PSA$Stage) <- c("pT2", "pT3")

#heatmap######################################################################
# heatmap data expression
Heat_data_boso <- boso_sujungtas_cleaned_PSA[, c("DN",
                                                 "PCA3_PSA_HPRT",
                                                 "PSMA_PSA_HPRT",
                                                 "AR_FL_PSA_HPRT"
                                                 
)]
Heat_data_boso <- as.data.frame(Heat_data_boso)
rownames(Heat_data_boso) <- Heat_data_boso[, 1]
Heat_data_boso <- Heat_data_boso[, -1]
colnames(Heat_data_boso) <- c("PCA3", "PSMA", "AR")
col_fun = colorRamp2(c(-5, -3, 0, 3, 5), c("#8564fb", "#64b3fb","white","#e088bd", "#af2745"))

heatmap <- Heatmap(as.matrix(Heat_data_boso), na_col = "grey", col = col_fun,
                   cluster_rows = FALSE,      # Disable row clustering
                   cluster_columns = FALSE)

plot(heatmap)

# heatmap data spd(l)s
Heat_data_boso_margaritos <- boso_sujungtas_cleaned_PSA[, c("DN",
                                                            "I sPDL1, pg/ml",
                                                            "I sPD1, pg/ml"
                                                            
)]
Heat_data_boso_margaritos <- as.data.frame(Heat_data_boso_margaritos)
rownames(Heat_data_boso_margaritos) <- Heat_data_boso_margaritos[, 1]
Heat_data_boso_margaritos <- Heat_data_boso_margaritos[, -1]
colnames(Heat_data_boso_margaritos) <- c("sPDL1", "sPD1")
col_fun_m  <- colorRampPalette(c("#e088bd", "#70a1fb"))(2)
heatmap_margaritos <- Heatmap(as.matrix(Heat_data_boso_margaritos), na_col = "grey",
                              col = col_fun_m,
                              cluster_rows = FALSE,      # Disable row clustering
                              cluster_columns = FALSE, 
                              column_names_gp = gpar(fontface = "italic"), 
                              name = "sPD(L) pg/mL", 
                              row_names_gp = gpar(fontsize = 8))

plot(heatmap_margaritos)

# clinical data as df
clinical_boso <- boso_sujungtas_cleaned_PSA[, c("DN", "clinical_significance" ,
                                                "grade_label"  ,"stage_number",
                                                "Amžius operacijos dienai",
                                                "PSA pries operacija ng/ml (Kokia tyrimo data???)" )]
colnames(clinical_boso) <- c("DN", "Clinical Significance", "Grade", "Stage", "Age", "PSA")
clinical_boso <- as.data.frame(clinical_boso)
rownames(clinical_boso) <- clinical_boso$patient_id
clinical_boso <- clinical_boso[, -1]
head(clinical_boso)

col_age <- colorRamp2(c(40, 90), c( "#9cd4c4", "#3c402f"))
col_PSA <- colorRamp2(c(0, 70), c( "#9cd4c4", "#3c402f"))
clinical_boso$Stage <- as.factor(clinical_boso$Stage)
levels(clinical_boso$Stage) <- c("pT2", "pT3")
# with clinical data

row_ha = rowAnnotation(Stage = clinical_boso$Stage, Grade = clinical_boso$Grade,
                       Age = clinical_boso$Age, PSA = clinical_boso$PSA,
                       col = list(Stage = c("pT2" = "#c8d49c",
                                            "pT3" = "#9cd4c4"), 
                                  Grade = c("grade 1" = "#9cd4c4",  
                                            "grade 2" = "#a89cd4", 
                                            "grade 3" = "#d49cac"), 
                                  Age = col_age,
                                  PSA = col_PSA
                       ))
labels <- c("\u22125", "0", "5") 

heatmap_raiska <- Heatmap(as.matrix(Heat_data_boso), na_col = "grey", col = col_fun,
                          cluster_rows = FALSE,      # Disable row clustering
                          cluster_columns = FALSE,right_annotation = row_ha,
                          row_split = clinical_boso$Stage,
                          column_names_gp = gpar(fontface = "italic"),
                          name = "Relative gene expression",heatmap_legend_param = list(
                            at = c(-5, 0, 5),   # Legend positions
                            labels = labels     # Adjusted labels
                          ))
plot(heatmap_raiska)

# list of heatmap
htlist <- heatmap_margaritos + heatmap_raiska

htlist 

draw(htlist, row_km = 1, row_split = clinical_boso$Stage)

#save png
png("heatmap_output0114.png", width = 3000, height = 2300, res = 300) # width and height in pixels, resolution in dpi
draw(htlist, row_km = 1, row_split = clinical_boso$Stage) # Render the heatmap
dev.off() # Close the PNG device

