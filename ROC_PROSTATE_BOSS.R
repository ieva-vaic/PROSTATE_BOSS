#Prostate cancer project of gene expression and sPD(L) - ROC analysis

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


#ROC: separate biomarkers###################################
roc_results_1<- lapply(main_biomarkers, function(col) {
  roc(response = boso_sujungtas_cleaned_PSA$Stage, predictor = boso_sujungtas_cleaned_PSA[[col]])
})
names(roc_results_1) <- main_biomarkers
roc_results_1
#extract the aucs
auc_values1 <- sapply(roc_results_1, function(roc_obj) {auc(roc_obj)})
auc_values1 
#figure
par(pty = "s") #sets square
plot.roc(roc_results_1[["PCA3_PSA_HPRT"]], print.auc = F, col = "#f77189") #
lines(roc_results_1[["PSMA_PSA_HPRT"]], col = "#4e9dd5", lwd =2) #
lines(roc_results_1[["AR_FL_PSA_HPRT"]], col ="#4eae30", lwd =2) #
lines(roc_results_1[["I sPDL1, pg/ml"]], col ="purple", lwd =2, lty =5) #
lines(roc_results_1[["I sPD1, pg/ml"]], col ="yellow3", lwd =2, lty =5) #
# Add legend
legend("bottomright", legend = c( expression(italic("PCA3")), 
                                  expression(italic("PSMA")),
                                  expression(italic("AR")),
                                  expression(italic("sPD-L1")),
                                  expression(italic("sPD-1"))
),
col = c("#f77189", "#4e9dd5","#4eae30", "purple", "yellow3" ), lty = 1, cex = 0.8, lwd =2)
title( main = "ROC Curves of Biomarker Predictors \n for Stage 3 PCa \n ", line = 1)
#table 
coords_results1 <- lapply(roc_results_1, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
})
coords_results1
# Create a dataframe combining AUC values and coordinates results
results_1<- data.frame(
  Predictor = c("PCA3", "PSMA", "AR", "sPD-L1", "sPD-1" ),
  AUC = auc_values1,
  do.call(rbind, coords_results1) 
)
rownames(results_1) <- NULL
results_1
#call names in lithuanian
colnames(results_1) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė",
                         "tikslumas", "jautrumas", "specifiškumas", 
                         "tpd", "npd", "tta", "nta")
#table nice design
formatted_1<- formattable(results_1, list(
  Biožymuo = formatter("span", 
                       style = ~ style(font.style = "italic"))))

formattable(formatted_1)
#nice table other way
gt_table1 <- results_1 %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai"
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table1

#logistic models: gene expression + spd(l)#############################
##gene expression + spdl ##############################################
logistic.model_2 <- glm(Stage~
                          PCA3_PSA_HPRT + PSMA_PSA_HPRT  + AR_FL_PSA_HPRT + `I sPDL1, pg/ml` ,
                        data = boso_sujungtas_cleaned_PSA,
                        family = "binomial")
logistic.model_2
pred_all_spdl1 <- predict.glm(logistic.model_2, type='response') #tik 2 predictions yra
summary(logistic.model_2)
pred_all_spdl1
# only leave the cases in the model
margaritos_pred_data2 <- boso_sujungtas_cleaned_PSA[(rownames(boso_sujungtas_cleaned_PSA) %in% names(pred_all_spdl1)), ]
table(response=margaritos_pred_data2$Stage)
dev_roc2<-roc(response=margaritos_pred_data2$Stage, predictor=pred_all_spdl1, data=margaritos_pred_data2)
auc(dev_roc2)
par(pty="s")
plot.roc(dev_roc2, main = "ROC: PSMA, PCA3 , AR + spdl1: kliniskai reiksmingi vs ne", print.auc=TRUE, print.thres=TRUE) 
#table 
coords_results2 <- 
  coords(dev_roc2, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_results2

##gene expression + spd ########################################################
logistic.model_3 <- glm(Stage~
                          PCA3_PSA_HPRT + PSMA_PSA_HPRT  + AR_FL_PSA_HPRT + `I sPD1, pg/ml` ,
                        data = boso_sujungtas_cleaned_PSA,
                        family = "binomial")
logistic.model_3
pred_all_spdl3 <- predict.glm(logistic.model_3, type='response') #tik 2 predictions yra
summary(logistic.model_3)
pred_all_spdl3
# only leave the cases in the model
margaritos_pred_data3 <- boso_sujungtas_cleaned_PSA[(rownames(boso_sujungtas_cleaned_PSA) %in% names(pred_all_spdl3)), ]
table(margaritos_pred_data3$Stage)
dev_roc3<-roc(response=margaritos_pred_data3$Stage, predictor=pred_all_spdl3, data=margaritos_pred_data3)
auc(dev_roc3)
par(pty="s")
plot.roc(dev_roc3, main = "ROC: PSMA, PCA3 , AR + spdl1: kliniskai reiksmingi vs ne", print.auc=TRUE, print.thres=TRUE) 
#table 
coords_results3<- 
  coords(dev_roc3, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_results3

##only gene expression model ###########################################################
logistic.model_4 <- glm(Stage~
                          PCA3_PSA_HPRT + PSMA_PSA_HPRT  + AR_FL_PSA_HPRT,
                        data = boso_sujungtas_cleaned_PSA,
                        family = "binomial")
logistic.model_4
pred_all_spdl4 <- predict.glm(logistic.model_4, type='response') #tik 2 predictions yra
summary(logistic.model_4)
pred_all_spdl4
# only leave the cases in the model
margaritos_pred_data4 <- boso_sujungtas_cleaned_PSA[(rownames(boso_sujungtas_cleaned_PSA) %in% names(pred_all_spdl4)), ]
table(margaritos_pred_data4$Stage)
dev_roc4<-roc(response=margaritos_pred_data4$Stage, 
              predictor=pred_all_spdl4, data=margaritos_pred_data4)
auc(dev_roc4)
#plot
par(pty="s")
plot.roc(dev_roc4, main = "ROC: PSMA, PCA3 , AR + spdl1: kliniskai reiksmingi vs ne", print.auc=TRUE, print.thres=TRUE) 
#table 
coords_results4<- 
  coords(dev_roc4, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_results4

##plot: gene expression + spd #####################################################
par(pty = "s") #sets square
plot.roc(dev_roc4, print.auc = F, col = "#f77189") #be margaritos
lines(dev_roc2, col = "#4e9dd5", lwd =2) # spdl1
lines(dev_roc3, col ="#4eae30", lwd =2) # spd
# Add legend
pastel_colors <- c("#f77189", "#4e9dd5","#4eae30" )
legend("bottomright", legend=c(expression(italic("PCA3, PSMA, AR")),
                               expression(italic("PCA3, PSMA, AR, sPDL1")),
                               expression(italic("PCA3, PSMA, AR, sPD1"))),
       col=pastel_colors, lty= 1, cex = 0.8)
title( main = "ROC Curves of Biomarker Predictor Combinantions \n for Stage 3 PCa \n ", line = 0.6)

#dataframe of coords
results_roc1<- data.frame(
  Biožymuo = c("PCA3, PSMA, AR", "PCA3, PSMA, AR, sPDL1", "PCA3, PSMA, AR, sPD1"),
  `plotas po kreive` = c(dev_roc4$auc, dev_roc2$auc, dev_roc3$auc), 
  `slenkstinė vertė` = c(coords_results4$threshold, coords_results2$threshold, coords_results3$threshold),
  tikslumas = c(coords_results4$accuracy , coords_results2$accuracy , coords_results3$accuracy),
  jautrumas = c(coords_results4$sensitivity, coords_results2$sensitivity, coords_results3$sensitivity),
  specifiškumas = c(coords_results4$specificity, coords_results2$specificity, coords_results3$specificity),
  tpd  = c(coords_results4$precision, coords_results2$precision, coords_results3$precision),
  npd  = c(coords_results4$npv, coords_results2$npv, coords_results3$npv),
  tta  = c(coords_results4$tpr, coords_results2$tpr, coords_results3$tpr),
  nta  = c(coords_results4$fpr, coords_results2$fpr, coords_results3$fpr),
  check.names = FALSE
)
rownames(results_roc1) <- NULL
results_roc1
#make nice table
gt_table <- results_roc1 %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai",
    subtitle = "Genų raiška + sPD-L1 arba sPD-1"
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table

#logistic models: spd(l)s + genes one by one###################################
##spd(l)s + psma #############################################################
logistic.model_5 <- glm(Stage~
                          `I sPD1, pg/ml` + `I sPDL1, pg/ml`  + PSMA_PSA_HPRT,
                        data = boso_sujungtas_cleaned_PSA,
                        family = "binomial")
logistic.model_5
pred_all_spdl5 <- predict.glm(logistic.model_5, type='response') #tik 2 predictions yra
summary(logistic.model_5)
pred_all_spdl5
# only leave the cases in the model
margaritos_pred_data5 <- boso_sujungtas_cleaned_PSA[(rownames(boso_sujungtas_cleaned_PSA) %in% names(pred_all_spdl5)), ]
table(margaritos_pred_data5$Stage)
dev_roc5<-roc(response=margaritos_pred_data5$Stage, 
              predictor=pred_all_spdl5, data=margaritos_pred_data5)
auc(dev_roc5)
#plot
par(pty="s")
plot.roc(dev_roc5, main = "ROC: PSMA, PCA3 , AR + spdl1: kliniskai reiksmingi vs ne", print.auc=TRUE, print.thres=TRUE) 
#table 
coords_results5<- 
  coords(dev_roc5, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_results5

##spd(l)s + PCA3#####################################
logistic.model_6 <- glm(Stage~
                          `I sPD1, pg/ml` + `I sPDL1, pg/ml`  + PCA3_PSA_HPRT,
                        data = boso_sujungtas_cleaned_PSA,
                        family = "binomial")
logistic.model_6
pred_all_spdl6 <- predict.glm(logistic.model_6, type='response') #tik 2 predictions yra
summary(logistic.model_6)
pred_all_spdl6
# only leave the cases in the model
margaritos_pred_data6 <- boso_sujungtas_cleaned_PSA[(rownames(boso_sujungtas_cleaned_PSA) %in% names(pred_all_spdl6)), ]
table(margaritos_pred_data6$Stage)
dev_roc6<-roc(response=margaritos_pred_data6$Stage, 
              predictor=pred_all_spdl6, data=margaritos_pred_data6)
auc(dev_roc6)
#plot
par(pty="s")
plot.roc(dev_roc6, main = "ROC: PSMA, PCA3 , AR + spdl1: kliniskai reiksmingi vs ne", print.auc=TRUE, print.thres=TRUE) 
#table 
coords_results6<- 
  coords(dev_roc6, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_results6

##spd(l)s +  + AR#####################################
logistic.model_7 <- glm(Stage~
                          `I sPD1, pg/ml` + `I sPDL1, pg/ml`  + AR_FL_PSA_HPRT,
                        data = boso_sujungtas_cleaned_PSA,
                        family = "binomial")
logistic.model_7
pred_all_spdl7 <- predict.glm(logistic.model_7, type='response') #tik 2 predictions yra
summary(logistic.model_7)
pred_all_spdl7
# only leave the cases in the model
margaritos_pred_data7 <- boso_sujungtas_cleaned_PSA[(rownames(boso_sujungtas_cleaned_PSA)
                                                     %in% names(pred_all_spdl7)), ]
table(margaritos_pred_data7$Stage)
dev_roc7<-roc(response=margaritos_pred_data7$Stage, 
              predictor=pred_all_spdl7, data=margaritos_pred_data7)
auc(dev_roc7)
#plot
par(pty="s")
plot.roc(dev_roc7, main = "ROC: PSMA, PCA3 , AR + spdl1: kliniskai reiksmingi vs ne", print.auc=TRUE, print.thres=TRUE) 
#table 
coords_results7<- 
  coords(dev_roc7, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_results7

##only spd(l)s#####################################
logistic.model_8 <- glm(Stage~
                          `I sPD1, pg/ml` + `I sPDL1, pg/ml`,
                        data = boso_sujungtas_cleaned_PSA,
                        family = "binomial")
logistic.model_8
pred_all_spdl8 <- predict.glm(logistic.model_8, type='response') #tik 2 predictions yra
summary(logistic.model_8)
pred_all_spdl8
# only leave the cases in the model
margaritos_pred_data8 <- boso_sujungtas_cleaned_PSA[(rownames(boso_sujungtas_cleaned_PSA)
                                                     %in% names(pred_all_spdl8)), ]
table(margaritos_pred_data8$Stage)
dev_roc8<-roc(response=margaritos_pred_data8$Stage, 
              predictor=pred_all_spdl8, data=margaritos_pred_data8)
auc(dev_roc8)
#plot
par(pty="s")
plot.roc(dev_roc8, main = "ROC: PSMA, PCA3 , AR + spdl1: kliniskai reiksmingi vs ne", print.auc=TRUE, print.thres=TRUE) 
#table 
coords_results8<- 
  coords(dev_roc8, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_results8

##plot: spd(l)s + gene expression one by one #####################################################
par(pty = "s") #sets square
plot.roc(dev_roc8, print.auc = F, col = "#f77189") #be margaritos
lines(dev_roc5, col = "#4e9dd5", lwd =2) # psma
lines(dev_roc6, col ="#4eae30", lwd =2) # pca3
lines(dev_roc7, col ="yellow3", lwd =2) # AR
# Add legend
pastel_colors <- c("#f77189", "#4e9dd5","#4eae30", "yellow3" )
legend("bottomright", legend=c(expression(italic("sPD-L1, sPD-1 ")),
                               expression(italic("PSMA, sPD-L1, sPD-1 ")),
                               expression(italic("PCA3, sPD-L1, sPD-1 ")),
                               expression(italic("AR, sPD-L1, sPD-1 "))),
       col=pastel_colors, lty= 1, cex = 0.8)
title( main = "ROC Curves of Biomarker Predictor Combinantions \n for Stage 3 PCa \n ", line = 0.6)

#dataframe of coords
results_roc2<- data.frame(
  Biožymuo = c("sPD-L1, sPD-1", "PSMA, sPD-L1, sPD-1 ","PCA3, sPD-L1, sPD-1 ", "AR, sPD-L1, sPD-1 " ),
  `plotas po kreive` = c(dev_roc8$auc, dev_roc5$auc, dev_roc6$auc, dev_roc7$auc), 
  `slenkstinė vertė` = c(coords_results8$threshold, coords_results5$threshold, coords_results6$threshold, coords_results7$threshold),
  tikslumas = c(coords_results8$accuracy , coords_results5$accuracy , coords_results6$accuracy, coords_results7$accuracy),
  jautrumas = c(coords_results8$sensitivity, coords_results5$sensitivity, coords_results6$sensitivity, coords_results7$sensitivity),
  spefiškumas = c(coords_results8$specificity, coords_results5$specificity, coords_results6$specificity, coords_results7$specificity),
  tpd  = c(coords_results8$precision, coords_results5$precision, coords_results6$precision, coords_results7$precision),
  npd  = c(coords_results8$npv, coords_results5$npv, coords_results6$npv,  coords_results7$npv),
  tta  = c(coords_results8$tpr, coords_results5$tpr, coords_results6$tpr, coords_results7$tpr),
  nta  = c(coords_results8$fpr, coords_results5$fpr, coords_results6$fpr, coords_results7$fpr)
  ,check.names = FALSE
)
rownames(results_roc2) <- NULL
results_roc2
#make nicer table
gt_table2 <- results_roc2 %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai",
    subtitle = "sPD-L1 ir sPD-1 + genų raiškos biožymenys"
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table2

#logistic models: AR gene expression + spd(l)s ########################################
## spd1 + ar#################################################
logistic.model_9 <- glm(Stage~
                          `I sPD1, pg/ml` +  AR_FL_PSA_HPRT,
                        data = boso_sujungtas_cleaned_PSA,
                        family = "binomial")
logistic.model_9
pred_all_spdl9 <- predict.glm(logistic.model_9, type='response') #tik 2 predictions yra
summary(logistic.model_9)
pred_all_spdl9
# only leave the cases in the model
margaritos_pred_data9 <- boso_sujungtas_cleaned_PSA[(rownames(boso_sujungtas_cleaned_PSA)
                                                     %in% names(pred_all_spdl9)), ]
table(margaritos_pred_data9$Stage)
dev_roc9<-roc(response=margaritos_pred_data9$Stage, 
              predictor=pred_all_spdl9, data=margaritos_pred_data9)
auc(dev_roc9)
#plot
par(pty="s")
plot.roc(dev_roc9, main = "ROC: PSMA, PCA3 , AR + spdl1: kliniskai reiksmingi vs ne", print.auc=TRUE, print.thres=TRUE) 
#table 
coords_results9<- 
  coords(dev_roc9, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_results9

## spl1 + ar ##############################################################
logistic.model_9.1 <- glm(Stage~
                            `I sPDL1, pg/ml` +  AR_FL_PSA_HPRT,
                          data = boso_sujungtas_cleaned_PSA,
                          family = "binomial")
logistic.model_9.1
pred_all_spdl9.1 <- predict.glm(logistic.model_9.1, type='response') #tik 2 predictions yra
summary(logistic.model_9.1)
pred_all_spdl9.1
# only leave the cases in the model
margaritos_pred_data9.1 <- boso_sujungtas_cleaned_PSA[(rownames(boso_sujungtas_cleaned_PSA)
                                                       %in% names(pred_all_spdl9.1)), ]
table(margaritos_pred_data9.1$Stage)
dev_roc9.1<-roc(response=margaritos_pred_data9.1$Stage, 
                predictor=pred_all_spdl9.1, data=margaritos_pred_data9.1)
auc(dev_roc9.1)
#plot
par(pty="s")
plot.roc(dev_roc9.1, main = "ROC: PSMA, PCA3 , AR + spdl1: kliniskai reiksmingi vs ne", print.auc=TRUE, print.thres=TRUE) 
#table 
coords_results9.1<- 
  coords(dev_roc9.1, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_results9.1

##plot: AR + spd(l)s#######################################
par(pty = "s") #sets square
plot.roc(dev_roc9, print.auc = F, col = "#f77189") #spd
lines(dev_roc9.1, col = "#4e9dd5", lwd =2) # spdl
lines(dev_roc7, col ="#4eae30", lwd =2) # kombinacija

# Add legend
pastel_colors <- c("#f77189", "#4e9dd5","#4eae30" )
legend("bottomright", legend=c(expression(italic("AR, sPD-1 ")),
                               expression(italic("AR, sPD-L1 ")),
                               expression(italic("AR, sPD-L1, sPD-1 "))),
       col=pastel_colors, lty= 1, cex = 0.8)
title( main = "ROC Curves of Biomarker Predictor Combinantions \n for Stage 3 PCa \n ", line = 0.6)

#dataframe of coords
results_roc9<- data.frame(
  Biožymuo = c("AR, sPD-1", "AR, sPD-L1 ","AR, sPD-L1, sPD-1 " ),
  `plotas po kreive` = c(dev_roc9$auc, dev_roc9.1$auc, dev_roc7$auc), 
  `slenkstinė vertė` = c(coords_results9$threshold, coords_results9.1$threshold,  coords_results7$threshold),
  tikslumas = c(coords_results9$accuracy , coords_results9.1$accuracy ,  coords_results7$accuracy),
  jautrumas = c(coords_results9$sensitivity, coords_results9.1$sensitivity, coords_results7$sensitivity),
  specifiškumas = c(coords_results9$specificity, coords_results9.1$specificity,  coords_results7$specificity),
  tpd  = c(coords_results9$precision, coords_results9.1$precision,  coords_results7$precision),
  npd  = c(coords_results9$npv, coords_results9.1$npv,   coords_results7$npv),
  tta  = c(coords_results9$tpr, coords_results9.1$tpr, coords_results7$tpr),
  nta  = c(coords_results9$fpr, coords_results9.1$fpr,  coords_results7$fpr)
  ,check.names = FALSE
)
rownames(results_roc9) <- NULL
results_roc9
#make nicer table
gt_table9 <- results_roc9 %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai",
    subtitle = "sPD-L1 arba sPD-1 + AR raiška "
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table9

#make figure 4: all models together#######################################
par(pty = "s") #sets square
plot.roc(dev_roc8, print.auc = F, col = "#f77189") #be margaritos
lines(dev_roc5, col = "#4e9dd5", lwd =2) # psma
lines(dev_roc6, col ="#4eae30", lwd =2) # pca3
lines(dev_roc7, col ="yellow3", lwd =2) # AR
lines(dev_roc9, col ="pink", lwd =2, lty = 2) # AR + spd
lines(dev_roc9.1, col ="deeppink", lwd =2, lty = 2) # AR + spdl
# Add legend
pastel_colors <- c("#f77189", "#4e9dd5","#4eae30", "yellow3", "pink", "deeppink" )
legend("bottomright", legend=c(expression(italic("sPD-L1, sPD-1 ")),
                               expression(italic("PSMA, sPD-L1, sPD-1 ")),
                               expression(italic("PCA3, sPD-L1, sPD-1 ")),
                               expression(italic("AR, sPD-L1, sPD-1 ")),
                               expression(italic("AR, sPD-1 ")),
                               expression(italic("AR, sPD-L1 "))),
       col=pastel_colors, lty= 1, cex = 0.7)
title( main = "ROC Curves of Biomarker Predictor Combinantions \n for Stage 3 PCa \n ", line = 0.6)

#dataframe of coords
results_roc2.2<- data.frame(
  Biožymuo = c("sPD-L1, sPD-1", "PSMA, sPD-L1, sPD-1 ","PCA3, sPD-L1, sPD-1 ",
               "AR, sPD-L1, sPD-1 ", "AR, sPD-1 ",  "AR, sPD-L1 "),
  `plotas po kreive` = c(dev_roc8$auc, dev_roc5$auc, dev_roc6$auc, dev_roc7$auc, 
                         dev_roc9$auc, dev_roc9.1$auc), 
  `slenkstinė vertė` = c(coords_results8$threshold, coords_results5$threshold,
                         coords_results6$threshold, coords_results7$threshold,
                         coords_results9$threshold, coords_results9.1$threshold),
  tikslumas = c(coords_results8$accuracy , coords_results5$accuracy ,
                coords_results6$accuracy, coords_results7$accuracy,
                coords_results9$accuracy, coords_results9.1$accuracy),
  jautrumas = c(coords_results8$sensitivity, coords_results5$sensitivity,
                coords_results6$sensitivity, coords_results7$sensitivity,
                coords_results9$sensitivity, coords_results9.1$sensitivity),
  specifiškumas = c(coords_results8$specificity, coords_results5$specificity,
                    coords_results6$specificity, coords_results7$specificity,
                    coords_results9$specificity, coords_results9.1$specificity),
  tpd  = c(coords_results8$precision, coords_results5$precision,
           coords_results6$precision, coords_results7$precision,
           coords_results9$precision, coords_results9.1$precision),
  ndp  = c(coords_results8$npv, coords_results5$npv,
           coords_results6$npv,  coords_results7$npv,
           coords_results9$npv,  coords_results9.1$npv),
  tta  = c(coords_results8$tpr, coords_results5$tpr,
           coords_results6$tpr, coords_results7$tpr,
           coords_results9$tpr, coords_results9.1$tpr),
  nta  = c(coords_results8$fpr, coords_results5$fpr,
           coords_results6$fpr, coords_results7$fpr,
           coords_results9$fpr, coords_results9.1$fpr),
  check.names = FALSE
  
)
rownames(results_roc2.2) <- NULL
results_roc2.2
#make nicer table
gt_table2.2 <- results_roc2.2 %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai",
    subtitle = "sPD-L1 ir sPD-1 + genų raiškos biožymenys"
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table2.2

