#Prostate cancer project of gene expression and sPD(L) - ROC analysis 2
#clinically significant vs not

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
library(webshot2)
library(pagedown)
setwd("~/rprojects/OTHER DATA/BOSO data")

#load data################################################
boso_sujungtas_cleaned_PSA <- readRDS("boso_sujungtas_cleaned_PSA.RDS")
table(boso_sujungtas_cleaned_PSA$clinical_significance)
table(boso_sujungtas_cleaned_PSA$clinical_significance_f)
#ROC: separate biomarkers ############################################
PSA_hprt_values3 <- c("PCA3_PSA_HPRT", "PSMA_PSA_HPRT", "AR_FL_PSA_HPRT", "I sPDL1, pg/ml", "I sPD1, pg/ml" )
roc_results_1<- lapply(PSA_hprt_values3, function(col) {
  roc(response = boso_sujungtas_cleaned_PSA$clinical_significance, predictor = boso_sujungtas_cleaned_PSA[[col]])
})
names(roc_results_1) <- PSA_hprt_values3
roc_results_1
#extract the aucs
auc_values1 <- sapply(roc_results_1, function(roc_obj) {auc(roc_obj)})
auc_values1 
#plot
par(pty = "s") #sets square
plot.roc(roc_results_1[["PCA3_PSA_HPRT"]], print.auc = F, col = "#f77189") #
lines(roc_results_1[["PSMA_PSA_HPRT"]], col = "#4e9dd5", lwd =2) #
lines(roc_results_1[["AR_FL_PSA_HPRT"]], col ="#4eae30", lwd =2)
lines(roc_results_1[["I sPDL1, pg/ml"]], col ="purple", lwd =2)
lines(roc_results_1[["I sPD1, pg/ml"]], col ="orange", lwd =2)#
# Add legend
legend("bottomright", legend = c( expression(italic("PCA3")), 
                                  expression(italic("PSMA")),
                                  expression(italic("AR")),
                                  expression(italic("sPD-L1")),
                                  expression(italic("sPD-1"))),
       col = c("#f77189", "#4e9dd5","#4eae30", "purple", "orange" ), lty = 1, cex = 0.8, lwd =2)
title( main = "ROC Curves of Biomarker Predictors \n for Clinically Significant PCa \n ", line = 1)
#table 
coords_results1 <- lapply(roc_results_1, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
})
coords_results1
# Create a dataframe combining AUC values and coordinates results
results_1<- data.frame(
  Predictor = c("PCA3", "PSMA", "AR", "sPD-L1", "sPD-1"),
  AUC = auc_values1,
  do.call(rbind, coords_results1) 
)
rownames(results_1) <- NULL
results_1

#table nice
formatted_1<- formattable(results_1, list(
  Predictor = formatter("span", 
                        style = ~ style(font.style = "italic"))))

formattable(formatted_1)
#nice table other way
gt_table1 <- results_1 %>%
  gt() %>%
  tab_header(
    title = "ROC metrics",
    subtitle = "Single biomarker predictions of clinically significant PCa"
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table1

#logistic models: gene expression + spd(l) ###############################################
boso_sujungtas_cleaned_PSA$clinical_significance_f <- as.numeric(as.factor(boso_sujungtas_cleaned_PSA$clinical_significance))-1
## genes +spdl1############################################################################
logistic.model_2 <- glm(clinical_significance_f ~
                          PCA3_PSA_HPRT + PSMA_PSA_HPRT  + AR_FL_PSA_HPRT + `I sPDL1, pg/ml` ,
                        data = boso_sujungtas_cleaned_PSA,
                        family = "binomial")
logistic.model_2
pred_all_spdl1 <- predict.glm(logistic.model_2, type='response') #tik 2 predictions yra
summary(logistic.model_2)
pred_all_spdl1
# only leave the cases in the model
margaritos_pred_data2 <- boso_sujungtas_cleaned_PSA[(rownames(boso_sujungtas_cleaned_PSA) %in% names(pred_all_spdl1)), ]
table(margaritos_pred_data2$clinical_significance_f)
dev_roc2<-roc(response=margaritos_pred_data2$clinical_significance_f, predictor=pred_all_spdl1, data=margaritos_pred_data2)
auc(dev_roc2)
#plot
par(pty="s")
plot.roc(dev_roc2, main = "ROC: PSMA, PCA3 , AR + sPD-L1: kliniskai reiksmingi vs ne", print.auc=TRUE, print.thres=TRUE) 
#table 
coords_results1 <- 
  coords(dev_roc2, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_results1


## genes +spd1#######################################################################
logistic.model_1 <- glm(clinical_significance_f~
                          PCA3_PSA_HPRT + PSMA_PSA_HPRT  + AR_FL_PSA_HPRT + `I sPD1, pg/ml` ,
                        data = boso_sujungtas_cleaned_PSA,
                        family = "binomial")
logistic.model_1
pred_all_spdl11 <- predict.glm(logistic.model_1, type='response') #tik 2 predictions yra
summary(logistic.model_1)
pred_all_spdl11
# only leave the cases in the model
margaritos_pred_data11 <- boso_sujungtas_cleaned_PSA[(rownames(boso_sujungtas_cleaned_PSA) %in% names(pred_all_spdl11)), ]
table(margaritos_pred_data11$clinical_significance_f)
dev_roc11<-roc(response=margaritos_pred_data11$clinical_significance_f, predictor=pred_all_spdl11, data=margaritos_pred_data11)
auc(dev_roc11)
#plot
par(pty="s")
plot.roc(dev_roc11, main = "ROC: PSMA, PCA3 , AR + spd1: kliniskai reiksmingi vs ne", print.auc=TRUE, print.thres=TRUE) 
#table 
coords_results11 <- 
  coords(dev_roc11, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_results11

## gene biomarkers only####################################################
logistic.model_3 <- glm(clinical_significance_f~
                          PCA3_PSA_HPRT + PSMA_PSA_HPRT  + AR_FL_PSA_HPRT  ,
                        data = boso_sujungtas_cleaned_PSA,
                        family = "binomial")
logistic.model_3
pred_all_spdl3 <- predict.glm(logistic.model_3, type='response') #tik 2 predictions yra
summary(logistic.model_3)
pred_all_spdl3
# only leave the cases in the model
margaritos_pred_data3 <- boso_sujungtas_cleaned_PSA[(rownames(boso_sujungtas_cleaned_PSA) %in% names(pred_all_spdl3)), ]
table(margaritos_pred_data3$clinical_significance_f)
dev_roc3<-roc(response=margaritos_pred_data3$clinical_significance_f, predictor=pred_all_spdl3, data=margaritos_pred_data3)
auc(dev_roc3)
#plot
par(pty="s")
plot.roc(dev_roc3, main = "ROC: PSMA, PCA3 , AR: kliniskai reiksmingi vs ne", print.auc=TRUE, print.thres=TRUE) 
#table 
coords_results3 <- 
  coords(dev_roc3, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_results3

##plot: gene expression + spd(l)s##############################################
par(pty = "s") #sets square
plot.roc(dev_roc3, print.auc = F, col = "#f77189") #
lines(dev_roc2, col = "#4e9dd5", lwd =2) #
lines(dev_roc11, col ="#4eae30", lwd =2) #
# Add legend
pastel_colors <- c("#f77189", "#4e9dd5","#4eae30" )
legend("bottomright", legend=c(expression(italic("PCA3, PSMA, AR")),
                               expression(italic("PCA3, PSMA, AR, sPD-L1")),
                               expression(italic("PCA3, PSMA, AR, sPD-1"))),
       col=pastel_colors, lty= 1, cex = 0.8)
title( main = "ROC Curves of Biomarker Combinantions \n for Clinically Significant PCa \n ", line = 0.6)

#dataframe of coords, roc2
results_roc2<- data.frame(
  Predictor = c("PCA3, PSMA, AR", "PCA3, PSMA, AR, sPD-L1", "PCA3, PSMA, AR, sPD-1"),
  AUC = c(dev_roc3$auc, dev_roc2$auc, dev_roc11$auc), 
  threshold = c(coords_results3$threshold, coords_results1$threshold, coords_results11$threshold),
  accuracy = c(coords_results3$accuracy , coords_results1$accuracy , coords_results11$accuracy),
  sensitivity = c(coords_results3$sensitivity, coords_results1$sensitivity, coords_results11$sensitivity),
  specificity = c(coords_results3$specificity, coords_results1$specificity, coords_results11$specificity),
  precision  = c(coords_results3$precision, coords_results1$precision, coords_results11$precision),
  npv  = c(coords_results3$npv, coords_results1$npv, coords_results11$npv),
  tpr  = c(coords_results3$tpr, coords_results1$tpr, coords_results11$tpr),
  fpr  = c(coords_results3$fpr, coords_results1$fpr, coords_results11$fpr)
)
rownames(results_roc2) <- NULL
results_roc2
#nicer table
gt_table2 <- results_roc2 %>%
  gt() %>%
  tab_header(
    title = "ROC metrics",
    subtitle = "Gene expression + sPD-L1 or sPD-1"
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table2


#logistic models: spd(l)s + genes one by one###################################
#pca3 +spd1 + spdl1 #############################################################
logistic.model_x1 <- glm(clinical_significance_f~
                           PCA3_PSA_HPRT + `I sPD1, pg/ml` +`I sPDL1, pg/ml` ,
                         data = boso_sujungtas_cleaned_PSA,
                         family = "binomial")
logistic.model_x1
pred_all_x1 <- predict.glm(logistic.model_x1, type='response') #tik 2 predictions yra
summary(logistic.model_x1)
pred_all_x1
# only leave the cases in the model
margaritos_pred_data_x1 <- boso_sujungtas_cleaned_PSA[(rownames(boso_sujungtas_cleaned_PSA) %in% names(pred_all_x1)), ]
table(margaritos_pred_data_x1$clinical_significance_f)
dev_roc_x1<-roc(response=margaritos_pred_data_x1$clinical_significance_f, predictor=pred_all_x1, data=margaritos_pred_data_x1)
auc(dev_roc_x1)
#plot
par(pty="s")
plot.roc(dev_roc_x1, main = "ROC: PCA3 + spd1 + spdl1: kliniskai reiksmingi vs ne", print.auc=TRUE, print.thres=TRUE) 
dev_roc_x1
#table 
coords_results_x1 <- 
  coords(dev_roc_x1, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_results_x1

##psma +spd1 + spdl1################################################
logistic.model_x2 <- glm(clinical_significance_f~
                           PSMA_PSA_HPRT + `I sPD1, pg/ml` +`I sPDL1, pg/ml` ,
                         data = boso_sujungtas_cleaned_PSA,
                         family = "binomial")
logistic.model_x2
pred_all_x2 <- predict.glm(logistic.model_x2, type='response') #tik 2 predictions yra
summary(logistic.model_x2)
pred_all_x2
# only leave the cases in the model
margaritos_pred_data_x2 <- boso_sujungtas_cleaned_PSA[(rownames(boso_sujungtas_cleaned_PSA) %in% names(pred_all_x2)), ]
table(margaritos_pred_data_x2$clinical_significance_f)
dev_roc_x2<-roc(response=margaritos_pred_data_x2$clinical_significance_f, predictor=pred_all_x2, data=margaritos_pred_data_x2)
auc(dev_roc_x2)
#plot
par(pty="s")
plot.roc(dev_roc_x2, main = "ROC: PSMA + spd1 + spdl1: kliniskai reiksmingi vs ne", print.auc=TRUE, print.thres=TRUE) 
dev_roc_x2
#table 
coords_results_x2 <- 
  coords(dev_roc_x2, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_results_x2

##AR  +spd1 + spdl1###########################################
logistic.model_x3 <- glm(clinical_significance_f~
                           AR_FL_PSA_HPRT + `I sPD1, pg/ml` +`I sPDL1, pg/ml` ,
                         data = boso_sujungtas_cleaned_PSA,
                         family = "binomial")
logistic.model_x3
pred_all_x3 <- predict.glm(logistic.model_x3, type='response') #tik 2 predictions yra
summary(logistic.model_x3)
pred_all_x3
# only leave the cases in the model
margaritos_pred_data_x3 <- boso_sujungtas_cleaned_PSA[(rownames(boso_sujungtas_cleaned_PSA) %in% names(pred_all_x3)), ]
table(margaritos_pred_data_x3$clinical_significance_f)
dev_roc_x3<-roc(response=margaritos_pred_data_x3$clinical_significance_f, predictor=pred_all_x3, data=margaritos_pred_data_x3)
auc(dev_roc_x3)
#plot
par(pty="s")
plot.roc(dev_roc_x3, main = "ROC: AR + spd1 + spdl1: kliniskai reiksmingi vs ne", print.auc=TRUE, print.thres=TRUE) 
dev_roc_x3
#table 
coords_results_x3 <- 
  coords(dev_roc_x3, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_results_x3

##plot: spd(l)s + gene expression one by one ####################################
par(pty = "s") #sets square
plot.roc(dev_roc_x1, print.auc = F, col = "#f77189") #
lines(dev_roc_x2, col = "#4e9dd5", lwd =2) #
lines(dev_roc_x3, col ="#4eae30", lwd =2) #
# Add legend
pastel_colors <- c("#f77189", "#4e9dd5","#4eae30" )
legend("bottomright", legend=c(expression(italic("PCA3, sPD-L1, sPD-1 ")),
                               expression(italic("PSMA, sPD-L1, sPD-1 ")),
                               expression(italic("AR, sPD-L1, sPD-1 "))),
       col=pastel_colors, lty= 1, cex = 0.8)
title( main = "ROC Curves of Biomarker Combinantions \n for Clinically Significant PCa \n  ", line = 0.6)

results_roc3<- data.frame(
  Predictor = c("PCA3, sPD-1, sPD-L1", "PSMA, sPD-1, sPD-L1", "AR, sPD-1, sPD-L1"),
  AUC = c(dev_roc_x1$auc, dev_roc_x2$auc, dev_roc_x3$auc), 
  threshold = c(coords_results_x1$threshold, coords_results_x2$threshold, coords_results_x3$threshold),
  accuracy = c(coords_results_x1$accuracy , coords_results_x2$accuracy , coords_results_x3$accuracy),
  sensitivity = c(coords_results_x1$sensitivity, coords_results_x2$sensitivity, coords_results_x3$sensitivity),
  specificity = c(coords_results_x1$specificity, coords_results_x2$specificity, coords_results_x3$specificity),
  precision  = c(coords_results_x1$precision, coords_results_x2$precision, coords_results_x3$precision),
  npv  = c(coords_results_x1$npv, coords_results_x2$npv, coords_results_x3$npv),
  tpr  = c(coords_results_x1$tpr, coords_results_x2$tpr, coords_results_x3$tpr),
  fpr  = c(coords_results_x1$fpr, coords_results_x2$fpr, coords_results_x3$fpr)
)
rownames(results_roc3) <- NULL
results_roc3

gt_table3 <- results_roc3 %>%
  gt() %>%
  tab_header(
    title = "ROC metrics",
    subtitle = "Gene expression + sPD-1 + sPD-L1"
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table3



#spd1 + spdl1 only ######################################
logistic.model_m <- glm(clinical_significance_f~
                          `I sPD1, pg/ml` +`I sPDL1, pg/ml` ,
                        data = boso_sujungtas_cleaned_PSA,
                        family = "binomial")
logistic.model_m
pred_all_m <- predict.glm(logistic.model_m, type='response') #tik 2 predictions yra
summary(logistic.model_m)
pred_all_m
margaritos_pred_data_m <- boso_sujungtas_cleaned_PSA[(rownames(boso_sujungtas_cleaned_PSA) %in% names(pred_all_m)), ]
table(margaritos_pred_data_m$clinical_significance_f)
dev_roc_m<-roc(response=margaritos_pred_data_m$clinical_significance_f, predictor=pred_all_m, data=margaritos_pred_data_m)
auc(dev_roc_m)
#plot
par(pty="s")
plot.roc(dev_roc_m, main = "ROC: spd1 + spdl1: kliniskai reiksmingi vs ne", print.auc=TRUE, print.thres=TRUE) 
dev_roc_m
#table 
coords_results_m <- 
  coords(dev_roc_m, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_results_m

