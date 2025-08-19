#Prostate cancer project of gene expression and sPD(L) - boxplots 

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

#read rds
boso_sujungtas_cleaned_PSA <- readRDS("boso_sujungtas_cleaned_PSA.RDS")

#boxplot: clinicaly significant vs not ######################################
biomarkers <- c("PCA3_PSA_HPRT" ,"PSMA_PSA_HPRT" ,
                "AR_FL_PSA_HPRT" ,  "I sPDL1, pg/ml","I sPD1, pg/ml"  )
lapply(boso_sujungtas_cleaned_PSA[biomarkers], shapiro.test)  #none is normal!

clin_sig <- melt(boso_sujungtas_cleaned_PSA, id.vars="clinical_significance",  measure.vars=biomarkers) #be outliers

clin_sig <- clin_sig %>%
  # Rename 4 to 4wd, f to Front, r to Rear
  mutate(variable = recode(variable, "PCA3_PSA_HPRT" = "PCA3"
                           , "PSMA_PSA_HPRT" = "PSMA"
                           , "AR_FL_PSA_HPRT" = "AR",
                           "I sPDL1, pg/ml" = "sPD-L1",
                           "I sPD1, pg/ml" = "sPD-1")) 
stat.test_clin <- clin_sig %>%
  group_by(variable) %>%
  pairwise_wilcox_test(value ~ clinical_significance, 
                       ref.group = "Clinically insignificant" ,
                       p.adjust.method = "BH")
stat.test_clin

#plot
## p-value brackets with comparisons to a reference sample
each.vs.ref_sig <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "Clinically insignificant",   "Clinically significant",   0.217 , 6, "PCA3",
  "Clinically insignificant",   "Clinically significant",     0.039 	,  4, "PSMA", #keicias
  "Clinically insignificant",   "Clinically significant",    0.232 	, 4, "AR",
  "Clinically insignificant",   "Clinically significant",    0.033   	, 75, "sPD-L1",#keicias
  "Clinically insignificant",   "Clinically significant",    0.818  	,75, "sPD-1"
)

x1 <- ggplot(clin_sig, aes(x=clinical_significance, y=value, color=variable)) +
  geom_boxplot( outlier.shape = NA, size = 1) +
  geom_jitter(color="black", size=1, alpha=0.5) + 
  ylab(NULL)+
  facet_wrap(.~ variable , scales = "free", nrow = 2, strip.position = "left", 
             labeller = labeller(variable= as_labeller(c("AR" = "Relative gene expression",
                                                         "PCA3" = "Relative gene expression",
                                                         "PSMA" = "Relative gene expression",
                                                         "sPD-1" = "sPD-1, pg/mL", "sPD-L1" = "sPD-L1, pg/mL")))) +
  add_pvalue(each.vs.ref_sig)+
  theme_classic()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.text = element_text(face = "italic"), 
    strip.placement = "outside", 
    strip.background = element_blank(),
    panel.background = element_rect(fill = NA, colour = 'black'),
    panel.grid = element_blank())+
  scale_x_discrete(labels = c("Clinically insigificant PCa", "Clinically sigificant PCa")) +
  stat_boxplot(geom ='errorbar')+
  scale_color_manual(values = c("#f77189", "#4eae30", "#4e9dd5", "grey", "black"))+
  guides(color=guide_legend(title="Biomarkers"))+ labs(x = NULL)

x1

# Get gtable object
g <- ggplotGrob(x1)
# Sometimes helpful to get an idea about the grobs' position
grid.show.layout(gtable:::gtable_layout(g)) 
# Add an extra top row (make some space)
g <- gtable_add_rows(x = g, heights = unit(0.65, 'cm'), pos = 2)
g <- gtable_add_rows(x = g, heights = unit(1, 'cm'), pos = 11)
# First strip
g <- gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = NA)),
                                  textGrob(label = "AR", 
                                           gp = gpar(col = "black",
                                                     fontface = "italic"))),
                     t = 3, l = 9, b = 3, r = 9, 
                     name = c("strip-top-2-rectg", "strip-top-2-text"))
# Second strip
g <- gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = NA)),
                                  textGrob(label = "PCA3", 
                                           gp = gpar(col = "black",
                                                     fontface = "italic"))),
                     t = 3, l = 15, b = 3, r = 15, 
                     name = c("strip-top-2-rectg", "strip-top-2-text")) 
g <- gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = NA)),
                                  textGrob(label = "PSMA", 
                                           gp = gpar(col = "black",
                                                     fontface = "italic"))),
                     t = 3, l = 21, b = 3, r = 21,
                     name = c("strip-top-3-rectg", "strip-top-3-text")) 
g <- gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = NA)),
                                  textGrob(label = "sPD-1", 
                                           gp = gpar(col = "black",
                                                     fontface = "italic"))),
                     t = 12, l = 9, b = 12, r = 9,
                     name = c("strip-top-4-rectg", "strip-top-4-text")) 
g <- gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = NA)),
                                  textGrob(label = "sPD-L1", 
                                           gp = gpar(col = "black",
                                                     fontface = "italic"))),
                     t = 12, l = 15, b = 12, r = 15,
                     name = c("strip-top-5-rectg", "strip-top-5-text"))
grid.newpage()

grid.draw(g)



#boxplot: stage#####################################################
lapply(boso_sujungtas_cleaned_PSA[biomarkers], shapiro.test)  #none is normal!

stage_rez <- melt(boso_sujungtas_cleaned_PSA, id.vars="stage_number",  measure.vars=biomarkers) #be outliers
stage_rez$stage_number <- as.factor(stage_rez$stage_number)
stage_rez <- stage_rez %>%
  # Rename 4 to 4wd, f to Front, r to Rear
  mutate(variable = recode(variable, "PCA3_PSA_HPRT" = "PCA3"
                           , "PSMA_PSA_HPRT" = "PSMA"
                           , "AR_FL_PSA_HPRT" = "AR",
                           "I sPDL1, pg/ml" = "sPD-L1",
                           "I sPD1, pg/ml" = "sPD-1")) 
stat.test_stage <- stage_rez %>%
  group_by(variable) %>%
  pairwise_wilcox_test(value ~ stage_number, ref.group = "2" )
stat.test_stage
#chek:
wilcox.test(boso_sujungtas_cleaned_PSA$`I sPDL1, pg/ml` ~ boso_sujungtas_cleaned_PSA$Stadija_skaicius) #cia turbut israsyti kiti 
wilcox.test(boso_sujungtas_cleaned_PSA$`I sPDL1, pg/ml` ~ boso_sujungtas_cleaned_PSA$stage_number)
#plot
## p-value brackets with comparisons to a reference sample
each.vs.ref_sig2 <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "2",   "3",   0.031     , 6, "PCA3",
  "2",   "3",      0.0003	,  4, "PSMA", #new
  "2",   "3",     0.397         	, 4, "AR",
  "2",   "3",     0.031, 100, "sPD-L1", #new
  "2",   "3",     0.071	, 100, "sPD-1"
)

# iki 3 skaiciu po kablelio numazinimas
each.vs.ref_sig2$p.adj <- ifelse(each.vs.ref_sig2$p.adj < 0.001, 
                                 "<0.001", 
                                 paste0(sprintf("%.3f", as.numeric(each.vs.ref_sig2$p.adj))))

x2 <- ggplot(stage_rez, aes(x=stage_number, y=value, color=variable)) +
  geom_boxplot( outlier.shape = NA, size = 1) +
  geom_jitter(color="black", size=1, alpha=0.5) + 
  ylab(NULL)+
  facet_wrap(.~ variable , scales = "free", nrow = 2, strip.position = "left", 
             labeller = labeller(variable= as_labeller(c("AR" = "Relative gene expression",
                                                         "PCA3" = "Relative gene expression",
                                                         "PSMA" = "Relative gene expression",
                                                         "sPD-1" = "sPD-1, pg/mL", "sPD-L1" = "sPD-L1, pg/mL")))) +
  add_pvalue(each.vs.ref_sig2)+ #p val
  theme_classic()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.text = element_text(face = "italic"), 
    strip.placement = "outside", 
    strip.background = element_blank(),
    panel.background = element_rect(fill = NA, colour = 'black'),
    panel.grid = element_blank())+
  scale_x_discrete(labels = c("pT2", "pT3")) +
  stat_boxplot(geom ='errorbar')+
  scale_color_manual(values = c("#f77189", "#4eae30", "#4e9dd5", "grey", "black"))+
  guides(color=guide_legend(title="Biomarkers"))+ labs(x = NULL)

x2
# Get gtable object
g <- ggplotGrob(x2)
# Sometimes helpful to get an idea about the grobs' position
grid.show.layout(gtable:::gtable_layout(g)) 
# Add an extra top row (make some space)
g <- gtable_add_rows(x = g, heights = unit(0.65, 'cm'), pos = 2)
g <- gtable_add_rows(x = g, heights = unit(1, 'cm'), pos = 11)
# First strip
g <- gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = NA)),
                                  textGrob(label = "AR", 
                                           gp = gpar(col = "black",
                                                     fontface = "italic"))),
                     t = 3, l = 9, b = 3, r = 9, 
                     name = c("strip-top-2-rectg", "strip-top-2-text"))
# Second strip
g <- gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = NA)),
                                  textGrob(label = "PCA3", 
                                           gp = gpar(col = "black",
                                                     fontface = "italic"))),
                     t = 3, l = 15, b = 3, r = 15, 
                     name = c("strip-top-2-rectg", "strip-top-2-text")) 
g <- gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = NA)),
                                  textGrob(label = "PSMA", 
                                           gp = gpar(col = "black",
                                                     fontface = "italic"))),
                     t = 3, l = 21, b = 3, r = 21,
                     name = c("strip-top-3-rectg", "strip-top-3-text")) 
g <- gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = NA)),
                                  textGrob(label = "sPD-1", 
                                           gp = gpar(col = "black",
                                                     fontface = "italic"))),
                     t = 12, l = 9, b = 12, r = 9,
                     name = c("strip-top-4-rectg", "strip-top-4-text")) 
g <- gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = NA)),
                                  textGrob(label = "sPD-L1", 
                                           gp = gpar(col = "black",
                                                     fontface = "italic"))),
                     t = 12, l = 15, b = 12, r = 15,
                     name = c("strip-top-5-rectg", "strip-top-5-text"))
grid.newpage()

grid.draw(g)

#boxplot: grade######################################################
grade_rez <- melt(boso_sujungtas_cleaned_PSA, id.vars="grade_label",  measure.vars=biomarkers) #be outliers
grade_rez$grade_label <- as.factor(grade_rez$grade_label)
grade_rez <- grade_rez %>%
  # Rename 4 to 4wd, f to Front, r to Rear
  mutate(variable = recode(variable, "PCA3_PSA_HPRT" = "PCA3"
                           , "PSMA_PSA_HPRT" = "PSMA"
                           , "AR_FL_PSA_HPRT" = "AR",
                           "I sPDL1, pg/ml" = "sPD-L1",
                           "I sPD1, pg/ml" = "sPD-1")) 
#grade 1 vs 2
stat.test_grade<- grade_rez %>%
  filter(grade_label != "grade 3") %>%
  group_by(variable) %>%
  pairwise_wilcox_test(value ~ grade_label, ref.group = "grade 1" , p.adjust.method = "BH")
stat.test_grade #psma
#grade 2 vs 3
stat.test_grade<- grade_rez %>%
  filter(grade_label != "grade 1") %>%
  group_by(variable) %>%
  pairwise_wilcox_test(value ~ grade_label, ref.group = "grade 2" , p.adjust.method = "BH")
stat.test_grade #spdl1
#grade 3 vs 1
stat.test_grade<- grade_rez %>%
  filter(grade_label != "grade 2") %>%
  group_by(variable) %>%
  pairwise_wilcox_test(value ~ grade_label, ref.group = "grade 1" , p.adjust.method = "BH")
stat.test_grade #psma

## p-value brackets with comparisons to a reference sample
each.vs.ref_sig3 <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "grade 1",   "grade 2",   0.011, 4, "PSMA", #new
  "grade 2",   "grade 3",      0.026	,  100, "sPD-L1",#new
  "grade 1",   "grade 3",     0.056	, 6, "PCA3",
  "grade 1",   "grade 3",     0.005 , 6, "PSMA"
)

x3 <- ggplot(grade_rez, aes(x=grade_label, y=value, color=variable)) +
  geom_boxplot( outlier.shape = NA, size = 1) +
  geom_jitter(color="black", size=1, alpha=0.5) + 
  ylab(NULL)+
  facet_wrap(.~ variable , scales = "free", nrow = 2, strip.position = "left", 
             labeller = labeller(variable= as_labeller(c("AR" = "Relative gene expression",
                                                         "PCA3" = "Relative gene expression",
                                                         "PSMA" = "Relative gene expression",
                                                         "sPD-1" = "sPD-1, pg/mL", "sPD-L1" = "sPD-L1, pg/mL")))) +
  add_pvalue(each.vs.ref_sig3)+ #p val
  theme_classic()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.text = element_text(face = "italic"), 
    strip.placement = "outside", 
    strip.background = element_blank(),
    panel.background = element_rect(fill = NA, colour = 'black'),
    panel.grid = element_blank())+
  scale_x_discrete(labels = c("Grade 1", "Grade 2", "Grade 3")) +
  stat_boxplot(geom ='errorbar')+
  scale_color_manual(values = c("#f77189", "#4eae30", "#4e9dd5", "grey", "black"))+
  guides(color=guide_legend(title="Biomarkers"))+ labs(x = NULL)

x3
# Get gtable object
g <- ggplotGrob(x3)
# Sometimes helpful to get an idea about the grobs' position
grid.show.layout(gtable:::gtable_layout(g)) 
# Add an extra top row (make some space)
g <- gtable_add_rows(x = g, heights = unit(0.65, 'cm'), pos = 2)
g <- gtable_add_rows(x = g, heights = unit(1, 'cm'), pos = 11)
# First strip
g <- gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = NA)),
                                  textGrob(label = "AR", 
                                           gp = gpar(col = "black",
                                                     fontface = "italic"))),
                     t = 3, l = 9, b = 3, r = 9, 
                     name = c("strip-top-2-rectg", "strip-top-2-text"))
# Second strip
g <- gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = NA)),
                                  textGrob(label = "PCA3", 
                                           gp = gpar(col = "black",
                                                     fontface = "italic"))),
                     t = 3, l = 15, b = 3, r = 15, 
                     name = c("strip-top-2-rectg", "strip-top-2-text")) 
g <- gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = NA)),
                                  textGrob(label = "PSMA", 
                                           gp = gpar(col = "black",
                                                     fontface = "italic"))),
                     t = 3, l = 21, b = 3, r = 21,
                     name = c("strip-top-3-rectg", "strip-top-3-text")) 
g <- gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = NA)),
                                  textGrob(label = "sPD-1", 
                                           gp = gpar(col = "black",
                                                     fontface = "italic"))),
                     t = 12, l = 9, b = 12, r = 9,
                     name = c("strip-top-4-rectg", "strip-top-4-text")) 
g <- gtable_add_grob(x = g,
                     grobs = list(rectGrob(gp = gpar(col = NA, 
                                                     fill = NA)),
                                  textGrob(label = "sPD-L1", 
                                           gp = gpar(col = "black",
                                                     fontface = "italic"))),
                     t = 12, l = 15, b = 12, r = 15,
                     name = c("strip-top-5-rectg", "strip-top-5-text"))
grid.newpage()

grid.draw(g)
