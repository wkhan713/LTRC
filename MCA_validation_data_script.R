---
title: "CH_project_03/08/2022"
author: "Wasay Khan"
date: "3/8/2022"
output: pdf_document
editor_options:
  chunk_output_type: console
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

### IMPORT DATA
```{r}

library(dplyr)

###------------------------NEW TWIST ASSAY DATA (07/2022)-----------------------###

twist <- readRDS("~/mca_project/twist_assay_data_val/combine_data_0721.rds")
twist$has_mca <- as.integer(twist$has_mca)
twist$has_snv <- as.integer(twist$has_snv)
twist <- twist[ ,c(1,2,3,4,5,6,7,8,9,10,11,22,42,82,83,84,14,86,23:41)]
twist <- twist %>% mutate(has_CH = case_when(has_mca == 0 ~ 0, has_snv == 0 ~ 0, has_mca == 1 ~ 1, has_snv == 1 ~ 1))
colnames(twist)[12] <- "has_mca_twist"
colnames(twist)[13] <- "has_snv_twist"
colnames(twist)[19] <- "has_CH_twist"
colnames(twist)[3] <- "Race"
colnames(twist)[4] <- "Ethnicity"
colnames(twist)[5] <- "birthDate"
colnames(twist)[6] <- "isDead"
colnames(twist)[7] <- "deathDate"
colnames(twist)[8] <- "lastRecordDate"
colnames(twist)[9] <- "DNADate"
colnames(twist)[10] <- "YAgeAtDNA"
colnames(twist)[11] <- "YAgeAtLastRecORDeath"
twist <- unique(twist)
#---#
#MCA UPDATE (Tara's data)
tara_mca <- read.table("~/mca_project/new_data/tara_Autosomal_MCA_to_BC_project_all_chrom_mca_calls.txt", header = TRUE)
tara_mca <- select(tara_mca,sample_id,Autosomal_MCA)
tara_mca <- unique(tara_mca)
tara_mca <- tara_mca %>% filter(Autosomal_MCA ==1)
colnames(tara_mca)[1] <- "GRID"
#---#
twist <- twist %>% left_join(tara_mca)
twist$has_mca_twist <- NULL
colnames(twist)[38] <- "has_mca_twist"
twist[, 38][is.na(twist[, 38])] <- 0
twist$has_CH_twist <- NULL
twist <- twist %>% mutate(has_CH = case_when(has_mca_twist == 0 ~ 0, has_snv_twist == 0 ~ 0, has_mca_twist == 1 ~ 1, has_snv_twist == 1 ~ 1))
#df <- select(twist,GRID,has_mca_twist,has_snv_twist,has_CH)
#count <- df %>% filter(has_mca_twist == 0 & has_snv_twist == 1 & has_CH == 0)

#remove sex with NA
twist <- twist[!is.na(twist$Sex_update),]
twist <- twist[!is.na(twist$Race),]
twist <- twist[!is.na(twist$Ethnicity),]
twist <- twist[!is.na(twist$Diabetes),]
twist <- twist[!is.na(twist$BMI_30),]

#turn snv columns to true and false
library(magrittr)
twist %<>% mutate_if(is.logical,as.numeric)






#ALL_CAUSE
data_all_cause <- read.table("~/mca_project/new_data/data_all_cause", header = TRUE)
data_all_cause_twist <- merge(twist,data_all_cause, by="GRID")
data_all_cause_twist <- data_all_cause_twist[ -c(38:47,50,51,54:56)]


#Additional step
#test <- select(data_all_cause_2,GRID,all_cause_cancer,all.cause_dna_to_diag,dna_to_lastRecord,has_mca,has_snv,CH)
test = data_all_cause_twist %>%
  mutate(Group = case_when(has_mca_twist ==1 & has_snv_twist ==0 ~ "has_mca",
                          has_snv_twist ==1 & has_mca_twist ==0 ~ "has_snv",
                          has_CH==1 ~ "ch",
                          TRUE ~ "ch-"))
test$all.cause_dna_to_diag[test$all_cause_cancer==0] = 15



```

```{r, fig.width=2, fig.height=2}

library(cmprsk)
combined_inc <- cuminc(ftime   = test$all.cause_dna_to_diag,  # failure time variable
                         fstatus = test$all_cause_cancer,  # variable with distinct codes for different causes of failure
                         group   = test$Group,  # estimates will calculated within groups
                         rho     = 0, # Power of the weight function used in the tests.
                         cencode = 0, # value of fstatus variable which indicates the failure time is censored.
                         )
## CIF and the variance of point estimates
timepoints(combined_inc, times = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13))

plot(combined_inc,col = 1:4,lwd = 2.5, lty = 1,ylab = "Cumulative Incidence of Hematologic Malignancies")
legend("topright",
       legend = c("CH+","MCA","SNV/indel","CH-"),
       cex=0.62,col = c("red","green","blue","black"),
       lty = 1,
       lwd = 2,
       bg = "white")



```


```{r}


#HR
library(survival)
library(survminer)
library(ggfortify)
test$PPM1D <- factor(test$PPM1D, levels = c(0,1))
#test$Group = relevel(test$Group, ref = "0")
hr <- coxph(Surv(all.cause_dna_to_diag, all_cause_cancer)~PPM1D, data  = test)
summary(hr)
exp(confint(hr))

cuminc_hr <- data.frame(all_cause=c('ASXL2','BRCC3','IDH2','SETBP1','CBL','JAK2','SF3B1','DNMT3A','SRSF2','GNAS','MPL','TET2','GNB1','NRAS','TP53','PPM1D'),
                 index=1:3,
                 HR=c(3.71,4.86,18.88,17.58,2.32,3.43,3.86,0.40,4.17,0.51,7.69,0.88,2.33,7.15,1.46),
                 lower=c(0.52,0.68,7.54,2.43,0.58,2.36,1.90,0.24,2.20,0.07,1.89,0.50,0.74,4.20,),
                 upper=c(26.50,34.80,47.29,127.01,9.35,4.98,7.84,7.88,3.67,313.25,1.53,7.27,12.18))
#view data
head(cuminc_hr)
#load ggplot2
library(ggplot2)
#create forest plot
ggplot(data=cuminc_hr, aes(y=index, x=HR, xmin=lower, xmax=upper, color=c("black","red","blue"))) +
  geom_point() +
  #xlim(c(-0.1,3)) +
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(cuminc_hr), labels=cuminc_hr$all_cause) +
  labs(x='HR (95% CI) Combined Malignancies') +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  theme_minimal()


#NEW AESTHETIC FOR HR
library(forestplot)
# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta <-
  structure(list(
    mean  = c(NA, NA, 5.28, 2.53, 8.31, 1),
    lower = c(NA, NA, 2.29, 1.85, 5.66, 0.999999),
    upper = c(NA, NA, 12.16, 3.47, 12.22, 1.000001)),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -6L),
    class = "data.frame")

#count <- test %>% filter(Group == "ch-")

tabletext<-cbind(
  c("", "", "MCA", "SNV/INDEL","CH+","CH- (reference)"),
  c("","Total","68","292","58","371"),
  c("","Event","49","133","43","68"),
  c("","Hazard Ratio \n 95% CI","5.28 (2.29-12.16)","2.53 (1.85-3.47)","8.31 (5.66-12.22)", NA),
  c("", "P-Value", "9.58e-05", "7.58e-09","<2.0e-16", NA))

forestplot(tabletext,
           cochrane_from_rmeta,
           graph.pos = 4,
           new_page = TRUE,
           is.summary=c(rep(TRUE, 2), rep(FALSE, 8), TRUE),
           clip=c(0,15.00),
           hrzl_lines = gpar(col = "#444444"),
           xlog=F,
           xlab="HR",
           txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex=0.65)),
           col=fpColors(box="black",
                        line="royalblue",
                        summary="red"))





























###-------------------------------- MERGE-------------------------------------- ###

#LYMPHOID
data_lymph <- read.table("~/mca_project/new_data/data_lymph", header = TRUE)
data_lymph_twist <- merge(twist,data_lymph, by="GRID")
data_lymph_twist <- data_lymph_twist[ -c(2,15,17:24,27:28,31:33)]

#MYELOID
data_myelo <- read.table("~/mca_project/new_data/data_myelo", header = TRUE)
data_myelo_twist <- merge(twist,data_myelo, by="GRID")
data_myelo_twist <- data_myelo_twist[ -c(2,15,17:24,27:28,31:33)]

#ALL_CAUSE
data_all_cause <- read.table("~/mca_project/new_data/data_all_cause", header = TRUE)
data_all_cause_twist <- merge(twist,data_all_cause, by="GRID")
data_all_cause_twist <- data_all_cause_twist[ -c(2,15,17:24,27:28,31:33)]

#CARDIOVASCULAR MORBIDITY
cardio_update <- read.table("~/mca_project/new_data/cardio_update", header = T)
cardio_twist <- merge(twist,cardio_update, by="GRID")
cardio_twist <- cardio_twist[ -c(2,15,17:24,46,47,55:58,62:66)]


###---------------------------- SURVIVAL ANALYSIS------------------------------ ###
library(survival)
library(survminer)
library(lubridate)
library(ggfortify)

#count <- data_lymph_twist %>% filter(lymphoid_cancer == 1 & has_mca_twist == 0 & has_snv_twist == 0)

#LYMPHOID
fit_lymphoid <- survfit(Surv(dna_to_lastRecord, lymphoid_cancer) ~ has_snv_twist + has_mca_twist, data = data_lymph_twist)
summary(fit_lymphoid)
CH.plot <- autoplot(fit_lymphoid, conf.int = FALSE)
CH.plot <- CH.plot +
  ggtitle("Lymphoid Malignancies Survival Analysis_TWIST") +
  labs(x = "Time(Years)", y = "Survival Probability") +
  guides(SCALE = "none") +
  #ylim(c(0, 1.0)) +
  labs(colour = "Legend") +
  scale_color_manual(labels = c("CH- | n=22","Only mCA | n=2,","Only SNV/indel | n=16,","CH+(mCA and SNV/indel) | n=5,"), values = c("red","dark green","blue","purple"))
print(CH.plot)

#MYELOID
fit_myeloid <- survfit(Surv(dna_to_lastRecord, myeloid_cancer) ~ has_snv_twist + has_mca_twist, data = data_myelo_twist)
summary(fit_myeloid)
CH.plot <- autoplot(fit_myeloid, conf.int = FALSE)
CH.plot <- CH.plot +
  ggtitle("Myeloid Malignancies Survival Analysis_TWIST") +
  labs(x = "Time(Years)", y = "Survival Probability") +
  guides(SCALE = "none") +
  #ylim(c(0, 1.0)) +
  labs(colour = "Legend") +
  scale_color_manual(labels = c("CH- | n=53","Only mCA | n=4,","Only SNV/indel | n=77,","CH+(mCA and SNV/indel) | n=39,"), values = c("red","dark green","blue","purple"))
print(CH.plot)

#ALL_CAUSE
fit_all_cause <- survfit(Surv(dna_to_lastRecord, all_cause_cancer) ~ has_snv_twist + has_mca_twist, data = data_all_cause_twist)
summary(fit_all_cause)
CH.plot <- autoplot(fit_all_cause, conf.int = FALSE)
CH.plot <- CH.plot +
  ggtitle("Combined Malignancies Survival Analysis_TWIST (Lymphoid+Myeloid") +
  labs(x = "Time(Years)", y = "Survival Probability") +
  guides(SCALE = "none") +
  #ylim(c(0, 1.0)) +
  labs(colour = "Legend") +
  scale_color_manual(labels = c("CH- | n=75","Only mCA | n=6,","Only SNV/indel | n=93,","CH+(mCA and SNV/indel) | n=44,"), values = c("red","dark green","blue","purple"))
print(CH.plot)

#CARDIO
fit_cardio_twist <- survfit(Surv(dna_to_lastRecord, cardio_mort_after_risk) ~ has_snv_twist + has_mca_twist, data = cardio_twist)
summary(fit_cardio_twist)
CH.plot <- autoplot(fit_cardio_twist, conf.int = FALSE)
CH.plot <- CH.plot +
  ggtitle("Cardiovascular Morbidity Survival Analysis_TWIST") +
  labs(x = "Time(Years)", y = "Survival Probability") +
  guides(SCALE = "none") +
  #ylim(c(0, 1.0)) +
  labs(colour = "Legend") +
  scale_color_manual(labels = c("CH- | n=292","Only mCA | n=9,","Only SNV/indel | n=242,","CH+(mCA and SNV/indel) | n=82,"), values = c("red","blue","black","brown"))
print(CH.plot)



```



### HAZARD  RATIO
```{r}

#ount <- data_lymph_twist %>% filter(lymphoid_cancer == 1 & has_mca_twist == 1 & has_snv_twist == 0)

#HR Forest Plot
#Fit a Cox proportional hazards model
hazard_lymph.cox <- coxph(Surv(as.numeric(dna_to_lastRecord), event=lymphoid_cancer) ~ has_mca_twist*has_snv_twist, data = data_lymph_twist)
summary(hazard_lymph.cox)
df_lymph_hr <- data.frame(Lymphoid=c('Only mCA (n=2)','Only SNV/indel (n=16)', 'CH+ (n=5)'),
                 index=1:3,
                 HR=c(8.34, 1.11, 1.94), # CH + confidence intervals might be wrong, ALL HAVE LOG VALUES! = ch-=0
                 lower=c(1.92, 0.58, 0.04),
                 upper=c(36.36, 2.11, 30.67))
#view data
head(df_lymph_hr)
#load ggplot2
library(ggplot2)
#create forest plot
ggplot(data=df_lymph_hr, aes(y=index, x=HR, xmin=lower, xmax=upper)) +
  geom_point() +
  #xlim(c(-0.1,3)) +
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(df_lymph_hr), labels=df_lymph_hr$Lymphoid) +
  labs(title='Lymphoid Malignancies', x='Hazard Ratio', y = '') +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  theme_minimal()


```

###Survival Analysis: MYELOID
```{r}

library(dplyr)


#count <- data_myelo_twist %>% filter(myeloid_cancer == 1 & has_mca_twist == 1 & has_snv_twist == 1)

#HR Forest Plot
#Fit a Cox proportional hazards model
hazard_myelo.cox <- coxph(Surv(as.numeric(dna_to_lastRecord), event=myeloid_cancer) ~ has_mca_twist*has_snv_twist, data = data_myelo_twist)
summary(hazard_myelo.cox)
df_myelo_hr <- data.frame(Myeloid=c('Only mCA (n=4)','Only SNV/indel (n=77)', 'CH+ (n=39)'),
                 index=1:3,
                 HR=c(9.33, 2.39, 9.36), # CH + confidence intervals might be wrong, ALL HAVE LOG VALUES! = ch-=0
                 lower=c(3.36, 1.69, 0.79),
                 upper=c(25.91, 3.40, 40.99))
#view data
head(df_myelo_hr)
#load ggplot2
library(ggplot2)
#create forest plot
ggplot(data=df_myelo_hr, aes(y=index, x=HR, xmin=lower, xmax=upper)) +
  geom_point() +
  #xlim(c(-0.1,3)) +
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(df_lymph_hr), labels=df_myelo_hr$Myeloid) +
  labs(title='Myeloid Malignancies', x='Hazard Ratio', y = '') +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  theme_minimal()


```

###Survival: ALL-CAUSE
```{r}


#count <- data_all_cause_twist %>% filter(all_cause_cancer == 1 & has_mca_twist == 1 & has_snv_twist == 0)

#HR Forest Plot
#Fit a Cox proportional hazards model
hazard_all_cause.cox <- coxph(Surv(as.numeric(dna_to_lastRecord), event=all_cause_cancer) ~ has_mca_twist*has_snv_twist, data = data_all_cause_twist)
summary(hazard_all_cause.cox)
df_all_cause_hr <- data.frame(all_cause=c('Only mCA (n=6)','Only SNV/indel (n=93)', 'CH+ (n=44)'),
                 index=1:3,
                 HR=c(8.19, 2.07, 7.12), # CH + confidence intervals might be wrong, ALL HAVE LOG VALUES! = ch-=0
                 lower=c(3.55, 1.53, 0.92),
                 upper=c(18.90, 2.81, 34.70))
#view data
head(df_all_cause_hr)
#load ggplot2
library(ggplot2)
#create forest plot
ggplot(data=df_all_cause_hr, aes(y=index, x=HR, xmin=lower, xmax=upper)) +
  geom_point() +
  #xlim(c(-0.1,3)) +
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(df_all_cause_hr), labels=df_all_cause_hr$all_cause) +
  labs(title='Combined Malignancies (Lymphoid + Myeloid)', x='Hazard Ratio', y = '') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_minimal()


```



##CARDIOVASCULAR MORTALITY (updated 03/28/2022): Import Cardiovascular data (Combined) All SNV
```{r}


#HAZARD RATIO
cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ has_snv_twist*has_mca_twist, data = cardio_twist)
summary(cardio_hr)

cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ smoke+has_snv_twist*has_mca_twist, data = cardio_twist)
summary(cardio_hr)

cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ bmi_30+has_snv_twist*has_mca_twist, data = cardio_twist)
summary(cardio_hr)

cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ age_30+has_snv_twist*has_mca_twist, data = cardio_twist)
summary(cardio_hr)

cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ Hyperlipidemia+has_snv_twist*has_mca_twist, data = cardio_twist)
summary(cardio_hr)

cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ Diabetes+has_snv_twist*has_mca_twist, data = cardio_twist)
summary(cardio_hr)

cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ Hypertension+has_snv_twist*has_mca_twist, data = cardio_twist)
summary(cardio_hr)



#count <- cardio_twist %>% filter(cardio_mort_after_risk == 1 & has_mca_twist == 1 & has_snv_twist == 0)

df_cardio_mort_hr <- data.frame(hr=c('CH+ (n=82)','Only SNV/Indel (n=242)','only mCA (n=9)','Smoking (n=422)','BMI>30 (n=170)','Age>30 (n=592)','Hyperlipidemia (n=365)','Diabetes (n=53)','Hypertension (n=80)'),
                 index=1:9,
                 HR=c(1.37,0.85,2.68,0.37,0.88,0.75,1.35,0.83,2.05),
                 lower=c(0.29,0.72,1.37,0.31,0.73,0.53,1.11,0.61,1.42),
                 upper=c(3.43,1.01,5.22,0.44,1.05,1.06,1.64,1.13,2.96))
#view data
head(df_cardio_mort_hr)
#load ggplot2
library(ggplot2)
#create forest plot
ggplot(data=df_cardio_mort_hr, aes(y=index, x=HR, xmin=lower, xmax=upper)) +
  geom_point() +
  #xlim(c(-0.1,3)) +
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(df_cardio_mort_hr), labels=df_cardio_mort_hr$hr) +
  labs(title='Cardiovascular Mortality', x='Hazard Ratio', y = '') +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  theme_minimal()











cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ smoke+bmi_30+age_30+Hyperlipidemia+Diabetes+Hypertension+has_snv_twist*has_mca_twist, data = cardio_twist)
summary(cardio_hr)


df_cardio_mort_hr <- data.frame(hr=c('CH+','Only SNV/Indel','only mCA','Smoking','BMI>30','Age>30','Hyperlipidemia','Diabetes','Hypertension'),
                 index=1:9,
                 HR=c(2.89,1.20,2.74,0.56,1.04,2.10,2.06,0.96,1.76),
                 lower=c(2.31,1.14,2.60,0.55,1.01,2.02,2.02,0.91,1.68),
                 upper=c(3.01,1.25,2.89,0.68,1.06,2.18,2.11,1.01,1.83))
#view data
head(df_cardio_mort_hr)
#load ggplot2
library(ggplot2)
#create forest plot
ggplot(data=df_cardio_mort_hr, aes(y=index, x=HR, xmin=lower, xmax=upper)) +
  geom_point() +
  #xlim(c(-0.1,3)) +
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(df_cardio_mort_hr), labels=df_cardio_mort_hr$hr) +
  labs(title='Cardiovascular Mortality', x='Hazard Ratio', y = '') +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  theme_minimal()




```




## HR
```{r}

cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ has_snv*has_mca, data = cardio_update)
summary(cardio_hr)

cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ has_snv*has_mca+mca_net, data = cardio_update)
summary(cardio_hr)

cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ smoke+has_snv*has_mca, data = cardio_update)
summary(cardio_hr)

cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ bmi_30+has_snv*has_mca, data = cardio_update)
summary(cardio_hr)

cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ age_30+has_snv*has_mca, data = cardio_update)
summary(cardio_hr)

cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ Hyperlipidemia+has_snv*has_mca, data = cardio_update)
summary(cardio_hr)

cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ Diabetes+has_snv*has_mca, data = cardio_update)
summary(cardio_hr)

cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ Hypertension+has_snv*has_mca, data = cardio_update)
summary(cardio_hr)

cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ JAK2_9PCNLOH+has_snv*has_mca, data = cardio_update)
summary(cardio_hr)

cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ DNMT3A+has_snv*has_mca, data = cardio_update)
summary(cardio_hr)

cardio_hr <- coxph(Surv(as.numeric(dna_to_lastRecord), event=cardio_mort_after_risk) ~ JAK2+has_snv*has_mca, data = cardio_update)
summary(cardio_hr)


df_cardio_mort_hr <- data.frame(hr=c('CH+ (n=188)','Only SNV/Indel (n=999)','DNMT3A (n=783)','JAK2 (n=226)','only mCA (n=2187)','CH+(JAK2_9PCNLOH) (n=103)','Smoking (n=24,846)','BMI>30 (n=11,934)','Age>30 (n=33,132)','Hyperlipidemia (n=13,965)','Diabetes (n=1,704)','Hypertension (n=3,490)'),
                 index=1:12,
                 HR=c(2.43,1.39,0.83,1.18,2.73,1.93,0.57,1.03,2.09,2.08,0.97,1.77),
                 lower=c(2.17,1.27,0.68,0.96,2.59,1.60,0.56,1.01,2.00,2.03,0.92,1.69),
                 upper=c(2.70,1.53,1.03,1.44,2.88,2.10,0.58,1.05,2.17,2.13,1.02,1.85))



#count <- cardio_update %>% filter(hassnv == 1, cardio_mort_after_risk == 1)


#view data
head(df_cardio_mort_hr)
#load ggplot2
library(ggplot2)
#create forest plot
ggplot(data=df_cardio_mort_hr, aes(y=index, x=HR, xmin=lower, xmax=upper)) +
  geom_point() +
  #xlim(c(-0.1,3)) +
  geom_errorbarh(height=.1) +
  scale_y_continuous(breaks=1:nrow(df_cardio_mort_hr), labels=df_cardio_mort_hr$hr) +
  labs(title='Cardiovascular Morbidity', x='Hazard Ratio', y = '') +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  theme_minimal()



```























# ALL SNV HR FOREST PLOT VALIDATION
```{r}


library(dplyr)

###------------------------NEW TWIST ASSAY DATA (07/2022)-----------------------###

twist <- readRDS("~/mca_project/twist_assay_data_val/combine_data_0721.rds")
twist$has_mca <- as.integer(twist$has_mca)
twist$has_snv <- as.integer(twist$has_snv)
twist <- twist[ ,c(1,2,3,4,5,6,7,8,9,10,11,22,42,82,83,84,14,86,23:41)]
twist <- twist %>% mutate(has_CH = case_when(has_mca == 0 ~ 0, has_snv == 0 ~ 0, has_mca == 1 ~ 1, has_snv == 1 ~ 1))
colnames(twist)[12] <- "has_mca_twist"
colnames(twist)[13] <- "has_snv_twist"
colnames(twist)[19] <- "has_CH_twist"
colnames(twist)[3] <- "Race"
colnames(twist)[4] <- "Ethnicity"
colnames(twist)[5] <- "birthDate"
colnames(twist)[6] <- "isDead"
colnames(twist)[7] <- "deathDate"
colnames(twist)[8] <- "lastRecordDate"
colnames(twist)[9] <- "DNADate"
colnames(twist)[10] <- "YAgeAtDNA"
colnames(twist)[11] <- "YAgeAtLastRecORDeath"
twist <- unique(twist)
#---#
#MCA UPDATE (Tara's data)
tara_mca <- read.table("~/mca_project/new_data/tara_Autosomal_MCA_to_BC_project_all_chrom_mca_calls.txt", header = TRUE)
tara_mca <- select(tara_mca,sample_id,Autosomal_MCA)
tara_mca <- unique(tara_mca)
tara_mca <- tara_mca %>% filter(Autosomal_MCA ==1)
colnames(tara_mca)[1] <- "GRID"
#---#
twist <- twist %>% left_join(tara_mca)
twist$has_mca_twist <- NULL
colnames(twist)[38] <- "has_mca_twist"
twist[, 38][is.na(twist[, 38])] <- 0
twist$has_CH_twist <- NULL
twist <- twist %>% mutate(has_CH = case_when(has_mca_twist == 0 ~ 0, has_snv_twist == 0 ~ 0, has_mca_twist == 1 ~ 1, has_snv_twist == 1 ~ 1))
#df <- select(twist,GRID,has_mca_twist,has_snv_twist,has_CH)
#count <- df %>% filter(has_mca_twist == 0 & has_snv_twist == 1 & has_CH == 0)

#remove sex with NA
twist <- twist[!is.na(twist$Sex_update),]
twist <- twist[!is.na(twist$Race),]
twist <- twist[!is.na(twist$Ethnicity),]
twist <- twist[!is.na(twist$Diabetes),]
twist <- twist[!is.na(twist$BMI_30),]

#turn snv columns to true and false
library(magrittr)
twist %<>% mutate_if(is.logical,as.numeric)


#ALL_CAUSE
data_all_cause <- read.table("~/mca_project/new_data/data_all_cause", header = TRUE)
data_all_cause_twist <- merge(twist,data_all_cause, by="GRID")
data_all_cause_twist <- data_all_cause_twist[ -c(38:47,50,51,54:56)]


#Additional step
#test <- select(data_all_cause_2,GRID,all_cause_cancer,all.cause_dna_to_diag,dna_to_lastRecord,has_mca,has_snv,CH)
test = data_all_cause_twist %>%
  mutate(Group = case_when(has_mca_twist ==1 & has_snv_twist ==0 ~ "has_mca",
                          has_snv_twist ==1 & has_mca_twist ==0 ~ "has_snv",
                          has_CH==1 ~ "ch",
                          TRUE ~ "ch-"))
test$all.cause_dna_to_diag[test$all_cause_cancer==0] = NA



```

```{r}

#select only ch- from gorup column
ref <- select(test,GRID,ASXL2,BRCC3,CBL,DNMT3A,GNAS,GNB1,IDH1,IDH2,JAK2,KRAS,MPL,NRAS,PPM1D,SETBP1,SF3B1,SRSF2,TET2,TP53,MPL,Group)
ref <-  ref %>% mutate(Group = replace(Group, Group == "has_snv", NA))
ref <-  ref %>% mutate(Group = replace(Group, Group == "has_mca", NA))
ref <-  ref %>% mutate(Group = replace(Group, Group == "ch", NA))
ref <-  ref %>% mutate(Group = replace(Group, Group == "ch-", 0))


###-------------------------------------------------------------------------SNV-----------------------------------------------------------------------------###

ref <-  ref %>% mutate(ASXL2 = replace(ASXL2, ASXL2 == 0, NA))
#merge 2 columns into 1
ref$ASXL2_hr <- paste(ref$ASXL2,ref$Group)
ref <-  ref %>% mutate(ASXL2_hr = replace(ASXL2_hr, ASXL2_hr == "NA NA", NA))
ref <-  ref %>% mutate(ASXL2_hr = replace(ASXL2_hr, ASXL2_hr == "NA 0", 0))
ref <-  ref %>% mutate(ASXL2_hr = replace(ASXL2_hr, ASXL2_hr == "1 NA", 1))
ref <- ref[ -c(2) ]


ref <-  ref %>% mutate(BRCC3 = replace(BRCC3, BRCC3 == 0, NA))
#merge 2 columns into 1
ref$BRCC3_hr <- paste(ref$BRCC3,ref$Group)
ref <-  ref %>% mutate(BRCC3_hr = replace(BRCC3_hr, BRCC3_hr == "NA NA", NA))
ref <-  ref %>% mutate(BRCC3_hr = replace(BRCC3_hr, BRCC3_hr == "NA 0", 0))
ref <-  ref %>% mutate(BRCC3_hr = replace(BRCC3_hr, BRCC3_hr == "1 NA", 1))
ref <- ref[ -c(2) ]

ref <-  ref %>% mutate(CBL = replace(CBL, CBL == 0, NA))
#merge 2 columns into 1
ref$CBL_hr <- paste(ref$CBL,ref$Group)
ref <-  ref %>% mutate(CBL_hr = replace(CBL_hr, CBL_hr == "NA NA", NA))
ref <-  ref %>% mutate(CBL_hr = replace(CBL_hr, CBL_hr == "NA 0", 0))
ref <-  ref %>% mutate(CBL_hr = replace(CBL_hr, CBL_hr == "1 NA", 1))
ref <- ref[ -c(2) ]

ref <-  ref %>% mutate(DNMT3A = replace(DNMT3A, DNMT3A == 0, NA))
#merge 2 columns into 1
ref$DNMT3A_hr <- paste(ref$DNMT3A,ref$Group)
ref <-  ref %>% mutate(DNMT3A_hr = replace(DNMT3A_hr, DNMT3A_hr == "NA NA", NA))
ref <-  ref %>% mutate(DNMT3A_hr = replace(DNMT3A_hr, DNMT3A_hr == "NA 0", 0))
ref <-  ref %>% mutate(DNMT3A_hr = replace(DNMT3A_hr, DNMT3A_hr == "1 NA", 1))
ref <- ref[ -c(2) ]

ref <-  ref %>% mutate(GNAS = replace(GNAS, GNAS == 0, NA))
#merge 2 columns into 1
ref$GNAS_hr <- paste(ref$GNAS,ref$Group)
ref <-  ref %>% mutate(GNAS_hr = replace(GNAS_hr, GNAS_hr == "NA NA", NA))
ref <-  ref %>% mutate(GNAS_hr = replace(GNAS_hr, GNAS_hr == "NA 0", 0))
ref <-  ref %>% mutate(GNAS_hr = replace(GNAS_hr, GNAS_hr == "1 NA", 1))
ref <- ref[ -c(2) ]

ref <-  ref %>% mutate(GNB1 = replace(GNB1, GNB1 == 0, NA))
#merge 2 columns into 1
ref$GNB1_hr <- paste(ref$GNB1,ref$Group)
ref <-  ref %>% mutate(GNB1_hr = replace(GNB1_hr, GNB1_hr == "NA NA", NA))
ref <-  ref %>% mutate(GNB1_hr = replace(GNB1_hr, GNB1_hr == "NA 0", 0))
ref <-  ref %>% mutate(GNB1_hr = replace(GNB1_hr, GNB1_hr == "1 NA", 1))
ref <- ref[ -c(2) ]

ref <-  ref %>% mutate(IDH1 = replace(IDH1, IDH1 == 0, NA))
#merge 2 columns into 1
ref$IDH1_hr <- paste(ref$IDH1,ref$Group)
ref <-  ref %>% mutate(IDH1_hr = replace(IDH1_hr, IDH1_hr == "NA NA", NA))
ref <-  ref %>% mutate(IDH1_hr = replace(IDH1_hr, IDH1_hr == "NA 0", 0))
ref <-  ref %>% mutate(IDH1_hr = replace(IDH1_hr, IDH1_hr == "1 NA", 1))
ref <- ref[ -c(2) ]

ref <-  ref %>% mutate(IDH2 = replace(IDH2, IDH2 == 0, NA))
#merge 2 columns into 1
ref$IDH2_hr <- paste(ref$IDH2,ref$Group)
ref <-  ref %>% mutate(IDH2_hr = replace(IDH2_hr, IDH2_hr == "NA NA", NA))
ref <-  ref %>% mutate(IDH2_hr = replace(IDH2_hr, IDH2_hr == "NA 0", 0))
ref <-  ref %>% mutate(IDH2_hr = replace(IDH2_hr, IDH2_hr == "1 NA", 1))
ref <- ref[ -c(2) ]

ref <-  ref %>% mutate(JAK2 = replace(JAK2, JAK2 == 0, NA))
#merge 2 columns into 1
ref$JAK2_hr <- paste(ref$JAK2,ref$Group)
ref <-  ref %>% mutate(JAK2_hr = replace(JAK2_hr, JAK2_hr == "NA NA", NA))
ref <-  ref %>% mutate(JAK2_hr = replace(JAK2_hr, JAK2_hr == "NA 0", 0))
ref <-  ref %>% mutate(JAK2_hr = replace(JAK2_hr, JAK2_hr == "1 NA", 1))
ref <- ref[ -c(2) ]

ref <-  ref %>% mutate(KRAS = replace(KRAS, KRAS == 0, NA))
#merge 2 columns into 1
ref$KRAS_hr <- paste(ref$KRAS,ref$Group)
ref <-  ref %>% mutate(KRAS_hr = replace(KRAS_hr, KRAS_hr == "NA NA", NA))
ref <-  ref %>% mutate(KRAS_hr = replace(KRAS_hr, KRAS_hr == "NA 0", 0))
ref <-  ref %>% mutate(KRAS_hr = replace(KRAS_hr, KRAS_hr == "1 NA", 1))
ref <- ref[ -c(2) ]

ref <-  ref %>% mutate(MPL = replace(MPL, MPL == 0, NA))
#merge 2 columns into 1
ref$MPL_hr <- paste(ref$MPL,ref$Group)
ref <-  ref %>% mutate(MPL_hr = replace(MPL_hr, MPL_hr == "NA NA", NA))
ref <-  ref %>% mutate(MPL_hr = replace(MPL_hr, MPL_hr == "NA 0", 0))
ref <-  ref %>% mutate(MPL_hr = replace(MPL_hr, MPL_hr == "1 NA", 1))
ref <- ref[ -c(2) ]

ref <-  ref %>% mutate(NRAS = replace(NRAS, NRAS == 0, NA))
#merge 2 columns into 1
ref$NRAS_hr <- paste(ref$NRAS,ref$Group)
ref <-  ref %>% mutate(NRAS_hr = replace(NRAS_hr, NRAS_hr == "NA NA", NA))
ref <-  ref %>% mutate(NRAS_hr = replace(NRAS_hr, NRAS_hr == "NA 0", 0))
ref <-  ref %>% mutate(NRAS_hr = replace(NRAS_hr, NRAS_hr == "1 NA", 1))
ref <- ref[ -c(2) ]

ref <-  ref %>% mutate(PPM1D = replace(PPM1D, PPM1D == 0, NA))
#merge 2 columns into 1
ref$PPM1D_hr <- paste(ref$PPM1D,ref$Group)
ref <-  ref %>% mutate(PPM1D_hr = replace(PPM1D_hr, PPM1D_hr == "NA NA", NA))
ref <-  ref %>% mutate(PPM1D_hr = replace(PPM1D_hr, PPM1D_hr == "NA 0", 0))
ref <-  ref %>% mutate(PPM1D_hr = replace(PPM1D_hr, PPM1D_hr == "1 NA", 1))
ref <- ref[ -c(2) ]

ref <-  ref %>% mutate(SETBP1 = replace(SETBP1, SETBP1 == 0, NA))
#merge 2 columns into 1
ref$SETBP1_hr <- paste(ref$SETBP1,ref$Group)
ref <-  ref %>% mutate(SETBP1_hr = replace(SETBP1_hr, SETBP1_hr == "NA NA", NA))
ref <-  ref %>% mutate(SETBP1_hr = replace(SETBP1_hr, SETBP1_hr == "NA 0", 0))
ref <-  ref %>% mutate(SETBP1_hr = replace(SETBP1_hr, SETBP1_hr == "1 NA", 1))
ref <- ref[ -c(2) ]

ref <-  ref %>% mutate(SF3B1 = replace(SF3B1, SF3B1 == 0, NA))
#merge 2 columns into 1
ref$SF3B1_hr <- paste(ref$SF3B1,ref$Group)
ref <-  ref %>% mutate(SF3B1_hr = replace(SF3B1_hr, SF3B1_hr == "NA NA", NA))
ref <-  ref %>% mutate(SF3B1_hr = replace(SF3B1_hr, SF3B1_hr == "NA 0", 0))
ref <-  ref %>% mutate(SF3B1_hr = replace(SF3B1_hr, SF3B1_hr == "1 NA", 1))
ref <- ref[ -c(2) ]

ref <-  ref %>% mutate(SRSF2 = replace(SRSF2, SRSF2 == 0, NA))
#merge 2 columns into 1
ref$SRSF2_hr <- paste(ref$SRSF2,ref$Group)
ref <-  ref %>% mutate(SRSF2_hr = replace(SRSF2_hr, SRSF2_hr == "NA NA", NA))
ref <-  ref %>% mutate(SRSF2_hr = replace(SRSF2_hr, SRSF2_hr == "NA 0", 0))
ref <-  ref %>% mutate(SRSF2_hr = replace(SRSF2_hr, SRSF2_hr == "1 NA", 1))
ref <- ref[ -c(2) ]

ref <-  ref %>% mutate(TET2 = replace(TET2, TET2 == 0, NA))
#merge 2 columns into 1
ref$TET2_hr <- paste(ref$TET2,ref$Group)
ref <-  ref %>% mutate(TET2_hr = replace(TET2_hr, TET2_hr == "NA NA", NA))
ref <-  ref %>% mutate(TET2_hr = replace(TET2_hr, TET2_hr == "NA 0", 0))
ref <-  ref %>% mutate(TET2_hr = replace(TET2_hr, TET2_hr == "1 NA", 1))
ref <- ref[ -c(2) ]

ref <-  ref %>% mutate(TP53 = replace(TP53, TP53 == 0, NA))
#merge 2 columns into 1
ref$TP53_hr <- paste(ref$TP53,ref$Group)
ref <-  ref %>% mutate(TP53_hr = replace(TP53_hr, TP53_hr == "NA NA", NA))
ref <-  ref %>% mutate(TP53_hr = replace(TP53_hr, TP53_hr == "NA 0", 0))
ref <-  ref %>% mutate(TP53_hr = replace(TP53_hr, TP53_hr == "1 NA", 1))
ref <- ref[ -c(2) ]



#combine with OG dataset
test <- cbind(test, ref)


```

```{r, fig.width=2, fig.height=2}

library(cmprsk)
combined_inc <- cuminc(ftime   = test$all.cause_dna_to_diag,  # failure time variable
                         fstatus = test$all_cause_cancer,  # variable with distinct codes for different causes of failure
                         group   = test$Group,  # estimates will calculated within groups
                         rho     = 0, # Power of the weight function used in the tests.
                         cencode = 0, # value of fstatus variable which indicates the failure time is censored.
                         )
## CIF and the variance of point estimates
timepoints(combined_inc, times = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13))

plot(combined_inc,col = 1:4,lwd = 2.5, lty = 1,ylab = "Cumulative Incidence of Hematologic Malignancies")
legend("topright",
       legend = c("CH+","MCA","SNV/indel","CH-"),
       cex=0.62,col = c("red","green","blue","black"),
       lty = 1,
       lwd = 2,
       bg = "white")



```


```{r}


#HR
library(survival)
library(survminer)
library(ggfortify)

test$DNMT3A <- factor(test$DNMT3A_hr, levels = c(0,1))
#test$Group = relevel(test$Group, ref = "0")
hr <- coxph(Surv(all.cause_dna_to_diag, all_cause_cancer)~PPM1D_hr, data  = test)
summary(hr)
exp(confint(hr))


#NEW AESTHETIC FOR HR
library(forestplot)
# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta <-
  structure(list(
    mean  = c(NA, NA, 1.17, 2.76, 23.39, 20.01, 20.97, 3.76, 6.51, 1.15, 1.73, 20.01, 8.22, 4.78, 0.78, 13.65, 1.67, 1 ),
    lower = c(NA, NA, 0.16, 0.37, 6.74, 2.33, 4.05, 2.39, 2.81, 0.67, 0.88, 2.33, 1.77, 2.50, 0.24, 6.30, 0.91, 0.9999999 ),
    upper = c(NA, NA, 8.52, 20.46, 81.20, 171.79, 108.53, 5.91, 15.06, 1.99, 3.38, 171.80, 38.11, 9.13, 2.50, 29.59, 3.05, 1.000001)),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -18L),
    class = "data.frame")

#df <- test[ -c(18:35)]
#count <- df %>% filter(MPL_hr == 0 & all_cause_cancer == 1)

tabletext<-cbind(
  c("", "","ASXL2","BRCC3","IDH2","SETBP1","CBL","JAK2","SF3B1","DNMT3A","SRSF2","GNAS","MPL","TET2","GNB1","TP53","PPM1D","CH-(reference)"),
  c("","Total","1","1","5","1","4","52","11","102","12","6","2","48","4","18","30","364"),
  c("","Event","1","1","5","1","2","33","8","16","10","1","2","13","3","15","13","68"),
  c("","Hazard Ratio \n 95% CI","1.17", "2.76", "23.39", "20.01", "20.97", "3.76", "6.51", "1.15", "1.73", "20.01", "8.22", "4.78", "0.78", "13.65", "1.67", NA),
  c("", "P-Value","0.88","0.32","6.9e-07","<0.05","2.86e-04","8.81e-09","1.22e-05","0.61","0.11","<0.05","<0.05","2.3e-06","0.68","3.57e-11","0.096",NA))

forestplot(tabletext,
           cochrane_from_rmeta,
           graph.pos = 6,
           new_page = TRUE,
           is.summary=c(rep(FALSE, 2)),
           clip=c(0,950.00),
           hrzl_lines = gpar(col = "#444444"),
           xlog=F,
           xlab="HR",
           txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex=0.70)),
           col=fpColors(box="black",
                        line="royalblue",
                        summary="blue"))




```





















# ALL MCA HR FOREST PLOT VALIDATION
```{r}


library(dplyr)

###------------------------NEW TWIST ASSAY DATA (07/2022)-----------------------###

twist <- readRDS("~/mca_project/twist_assay_data_val/combine_data_0721.rds")
twist$has_mca <- as.integer(twist$has_mca)
twist$has_snv <- as.integer(twist$has_snv)
twist <- twist[ ,c(1,2,3,4,5,6,7,8,9,10,11,22,42,82,83,84,14,86,23:41,112,109,103,94,143,108,126,95,123,97,133,98,120,124,105,114,121,129,92,99)]
twist <- twist %>% mutate(has_CH = case_when(has_mca == 0 ~ 0, has_snv == 0 ~ 0, has_mca == 1 ~ 1, has_snv == 1 ~ 1))
colnames(twist)[12] <- "has_mca_twist"
colnames(twist)[13] <- "has_snv_twist"
colnames(twist)[19] <- "has_CH_twist"
colnames(twist)[3] <- "Race"
colnames(twist)[4] <- "Ethnicity"
colnames(twist)[5] <- "birthDate"
colnames(twist)[6] <- "isDead"
colnames(twist)[7] <- "deathDate"
colnames(twist)[8] <- "lastRecordDate"
colnames(twist)[9] <- "DNADate"
colnames(twist)[10] <- "YAgeAtDNA"
colnames(twist)[11] <- "YAgeAtLastRecORDeath"
twist <- unique(twist)
#---#
#MCA UPDATE (Tara's data)
tara_mca <- read.table("~/mca_project/new_data/tara_Autosomal_MCA_to_BC_project_all_chrom_mca_calls.txt", header = TRUE)
tara_mca <- select(tara_mca,sample_id,Autosomal_MCA)
tara_mca <- unique(tara_mca)
tara_mca <- tara_mca %>% filter(Autosomal_MCA ==1)
colnames(tara_mca)[1] <- "GRID"
#---#
twist <- twist %>% left_join(tara_mca)
twist$has_mca_twist <- NULL
colnames(twist)[58] <- "has_mca_twist"
twist[, 58][is.na(twist[, 58])] <- 0
twist$has_CH_twist <- NULL
twist <- twist %>% mutate(has_CH = case_when(has_mca_twist == 0 ~ 0, has_snv_twist == 0 ~ 0, has_mca_twist == 1 ~ 1, has_snv_twist == 1 ~ 1))
#df <- select(twist,GRID,has_mca_twist,has_snv_twist,has_CH)
#count <- df %>% filter(has_mca_twist == 0 & has_snv_twist == 1 & has_CH == 0)

#remove sex with NA
twist <- twist[!is.na(twist$Sex_update),]
twist <- twist[!is.na(twist$Race),]
twist <- twist[!is.na(twist$Ethnicity),]
twist <- twist[!is.na(twist$Diabetes),]
twist <- twist[!is.na(twist$BMI_30),]

#turn snv columns to true and false
library(magrittr)
twist %<>% mutate_if(is.logical,as.numeric)


#ALL_CAUSE
data_all_cause <- read.table("~/mca_project/new_data/data_all_cause", header = TRUE)
data_all_cause_twist <- merge(twist,data_all_cause, by="GRID")
data_all_cause_twist <- data_all_cause_twist[ -c(59:67,70,71,74:76)]


#Additional step
#test <- select(data_all_cause_2,GRID,all_cause_cancer,all.cause_dna_to_diag,dna_to_lastRecord,has_mca,has_snv,CH)
test = data_all_cause_twist %>%
  mutate(Group = case_when(has_mca_twist ==1 & has_snv_twist ==0 ~ "has_mca",
                          has_snv_twist ==1 & has_mca_twist ==0 ~ "has_snv",
                          has_CH==1 ~ "ch",
                          TRUE ~ "ch-"))
test$all.cause_dna_to_diag[test$all_cause_cancer==0] = NA


```

```{r}

#select only ch- from gorup column
ref <- select(test,GRID,'9pCN-LOH','7qDEL','5qDEL','2pDEL','20qDEL','7pDEL','13qDEL','2qDEL','12pDEL','3pDEL','17pDEL','3qDEL','11pDEL','12qDEL','6pDEL','9qDEL','11qDEL','16pDEL','1pDEL','4qCN-LOH',Group)
ref <-  ref %>% mutate(Group = replace(Group, Group == "has_snv", NA))
ref <-  ref %>% mutate(Group = replace(Group, Group == "has_mca", NA))
ref <-  ref %>% mutate(Group = replace(Group, Group == "ch", NA))
ref <-  ref %>% mutate(Group = replace(Group, Group == "ch-", 0))


###-------------------------------------------------------------------------SNV-----------------------------------------------------------------------###

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "9pCN-LOH_hr"

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- '7qDEL_hr'

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- '5qDEL_hr'

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "2pDEL_hr"

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "20qDEL_hr"

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "7pDEL_hr"

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "13qDEL_hr"

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "2qDEL_hr"

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "12pDEL_hr"

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "3pDEL_hr"

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "17pDEL_hr"

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "3qDEL_hr"

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "11pDEL_hr"

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "12qDEL_hr"

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "6pDEL_hr"

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "9qDEL_hr"

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "11qDEL_hr"

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "16pDEL_hr"

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "1pDEL_hr"

colnames(ref)[2] <- "mca"
ref$mca <- gsub(0, NA, ref$mca)
#merge 2 columns into 1
ref$hr <- paste(ref$mca,ref$Group)
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA NA", NA))
ref <-  ref %>% mutate(hr = replace(hr, hr == "NA 0", 0))
ref <-  ref %>% mutate(hr = replace(hr, hr == "1 NA", 1))
ref <- ref[ -c(2) ]
colnames(ref)[22] <- "4qCN-LOH_hr"



#combine with OG dataset
test <- cbind(test, ref)


```

```{r, fig.width=2, fig.height=2}

library(cmprsk)
combined_inc <- cuminc(ftime   = test$all.cause_dna_to_diag,  # failure time variable
                         fstatus = test$all_cause_cancer,  # variable with distinct codes for different causes of failure
                         group   = test$Group,  # estimates will calculated within groups
                         rho     = 0, # Power of the weight function used in the tests.
                         cencode = 0, # value of fstatus variable which indicates the failure time is censored.
                         )
## CIF and the variance of point estimates
timepoints(combined_inc, times = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13))

plot(combined_inc,col = 1:4,lwd = 2.5, lty = 1,ylab = "Cumulative Incidence of Hematologic Malignancies")
legend("topright",
       legend = c("CH+","MCA","SNV/indel","CH-"),
       cex=0.62,col = c("red","green","blue","black"),
       lty = 1,
       lwd = 2,
       bg = "white")



```


```{r}


#HR
library(survival)
library(survminer)
library(ggfortify)

colnames(df)[66] <- 'mca_hr'
test$mca_hr <- factor(test$mca_hr, levels = c(0,1))
hr <- coxph(Surv(all.cause_dna_to_diag, all_cause_cancer)~ mca_hr, data  = test)
summary(hr)
exp(confint(hr))
#colnames(test)[66] <- 'done'


#NEW AESTHETIC FOR HR
library(forestplot)
# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta <-
  structure(list(
    mean  = c(NA, NA, 2.82,24.08,8.77,22.65,4.10,23.39,12.04,20.97,23.39,20.01,24.08,20.97,21.84,20.97,5.10,20.01,11.48,1),
    lower = c(NA, NA, 1.50,7.31,4.21,6.05,0.96,6.74,4.19,4.05,6.74,2.33,7.31,4.05,5.20,4.05,2.43,2.33,2.45,0.9999999),
    upper = c(NA, NA, 5.30,79.30,18.27,84.76,17.46,81.20,34.65,108.53,81.20,171.80,79.30,108.53,91.83,108.53,10.70,171.80,53.69,1.0000001)),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -20L),
    class = "data.frame")

#df <- test[ , c(1,59,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85)]
#colnames(df)[3] <- 'mca_hr'
#count <- df %>% filter(mca_hr == 1)
#count <- df %>% filter(mca_hr == 1 & all_cause_cancer == 1)
#df <- df[ , -c(3)]

tabletext<-cbind(
  c("", "MCA",'9pCN-LOH','7qDEL','5qDEL','2pDEL','20qDEL','7pDEL','13qDEL','2qDEL','12pDEL','17pDEL','11pDEL','12qDEL','6pDEL','9qDEL','11qDEL','1pDEL','4qCN-LOH','CH-(reference)'),
   c("","Total","15","7","12","4","2","5","5","2","5","1","6","2","3","2","9","1","3","364"),
   c("","Event","12","6","12","4","2","5","5","2","5","1","6","2","3","2","9","1","2","68"),
   c("","Hazard Ratio \n 95% CI",'2.82','24.08','8.77','22.65','4.10','23.39','12.04','20.97','23.39','20.01','24.08','20.97','21.84','20.97','5.10','20.01','11.48', NA),
  c("", "P-Value","0.001","1.68e-07","6.45e-09","3.6e-06","0.06","6.9e-07","3.93e-06","2.86e-04","6.9e-07","0.006","1.68e-07","2.86e-04","2.57e-05","2.86e-04","1.66e-05","0.006","0.002",NA))

#df <- test[ , c(1,59,61,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85)]
#colnames(df)[4] <- 'mca_hr'
#hr <- coxph(Surv(all.cause_dna_to_diag, all_cause_cancer)~ mca_hr, data  = df)
#summary(hr)
#df <- df[ , -c(4)]

forestplot(tabletext,
           cochrane_from_rmeta,
           graph.pos = 6,
           new_page = TRUE,
           is.summary=c(rep(FALSE, 2)),
           clip=c(0,175.00),
           hrzl_lines = gpar(col = "#444444"),
           xlog=F,
           xlab="HR",
           txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex=0.70)),
           col=fpColors(box="black",
                        line="royalblue",
                        summary="blue"))


```
