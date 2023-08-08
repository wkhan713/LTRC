---
title: "7965-DH:T-cell vs. temperature in GRCm39"
author: "Wasay Khan"
date: '2022-05-05'
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## 7965 Heintzman Data
```{r}

library(dplyr)

sample_1_4 <- read.table("~/7965/7965-DH-0001_0004.vcf")
colnames(sample_1_4) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "7965-DH-0001", "7965-DH-0004")
sample_1_4 <- filter(sample_1_4, FILTER == "PASS")
#Isolate values
library(magrittr)
#GENE
GENE<-data.frame(sample_1_4[[8]]%>%strsplit("|", fixed=T)%>%sapply("[", 4))
colnames(GENE)[1] <- "GENE"
#DP
DP<-data.frame(sample_1_4[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 7))
colnames(DP)[1] <- "DP"
#F1R2
F1R2<-data.frame(sample_1_4[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 5))
colnames(F1R2)[1] <- "F1R2"
#F2R1
F2R1 <-data.frame(sample_1_4[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 6))
colnames(F2R1)[1] <- "F2R1"
#genotype
GT <-data.frame(sample_1_4[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 1))
colnames(GT)[1] <- "GT"
#combine all isolted columns to original dataframe
sample_1_4 <- cbind(sample_1_4,GENE,DP,F1R2,F2R1,GT)
sample_1_4 <-sample_1_4[,c(1:7,12,8:11,13:16)]
#DP >50
sample_1_4 <-sample_1_4[sample_1_4$DP>50,]
#remove GT = 0|1
sample_1_4 <-sample_1_4[sample_1_4$GT != "0|1",]
#remove F1R2,F2R1 values that end it ",0" ensuring reads in both forward and reverse direction
sample_1_4 <- sample_1_4 %>% filter(!grepl(",0$", F1R2), !grepl(",0$", F2R1))
write.table(sample_1_4, file = "~/7965/7965-DH-0001_0004")

#sample_1_4 <- read.table("~/7965/7965-DH-0001_0004")



```

```{r}



sample_1_7 <- read.table("~/7965/7965-DH-0001_0007.vcf")
colnames(sample_1_7) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "7965-DH-0001", "7965-DH-0007")
sample_1_7 <- filter(sample_1_7, FILTER == "PASS")
#Isolate values
library(magrittr)
#GENE
GENE<-data.frame(sample_1_7[[8]]%>%strsplit("|", fixed=T)%>%sapply("[", 4))
colnames(GENE)[1] <- "GENE"
#DP
DP<-data.frame(sample_1_7[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 7))
colnames(DP)[1] <- "DP"
#F1R2
F1R2<-data.frame(sample_1_7[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 5))
colnames(F1R2)[1] <- "F1R2"
#F2R1
F2R1 <-data.frame(sample_1_7[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 6))
colnames(F2R1)[1] <- "F2R1"
#genotype
GT <-data.frame(sample_1_7[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 1))
colnames(GT)[1] <- "GT"
#combine all isolted columns to original dataframe
sample_1_7 <- cbind(sample_1_7,GENE,DP,F1R2,F2R1,GT)
sample_1_7 <-sample_1_7[,c(1:7,12,8:11,13:16)]
#DP >50
sample_1_7 <-sample_1_7[sample_1_7$DP>50,]
#remove GT = 0|1
sample_1_7 <-sample_1_7[sample_1_7$GT != "0|1",]
#remove F1R2,F2R1 values that end it ",0" ensuring reads in both forward and reverse direction
sample_1_7 <- sample_1_7 %>% filter(!grepl(",0$", F1R2), !grepl(",0$", F2R1))
write.table(sample_1_7, file = "~/7965/7965-DH-0001_0007")

```

```{r}


sample_2_5 <- read.table("~/7965/7965-DH-0002_0005.vcf")
colnames(sample_2_5) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "7965-DH-0002", "7965-DH-0005")
sample_2_5 <- filter(sample_2_5, FILTER == "PASS")
#Isolate values
library(magrittr)
#GENE
GENE<-data.frame(sample_2_5[[8]]%>%strsplit("|", fixed=T)%>%sapply("[", 4))
colnames(GENE)[1] <- "GENE"
#DP
DP<-data.frame(sample_2_5[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 7))
colnames(DP)[1] <- "DP"
#F1R2
F1R2<-data.frame(sample_2_5[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 5))
colnames(F1R2)[1] <- "F1R2"
#F2R1
F2R1 <-data.frame(sample_2_5[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 6))
colnames(F2R1)[1] <- "F2R1"
#genotype
GT <-data.frame(sample_2_5[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 1))
colnames(GT)[1] <- "GT"
#combine all isolted columns to original dataframe
sample_2_5 <- cbind(sample_2_5,GENE,DP,F1R2,F2R1,GT)
sample_2_5 <-sample_2_5[,c(1:7,12,8:11,13:16)]
#DP >50
sample_2_5 <-sample_2_5[sample_2_5$DP>50,]
#remove GT = 0|1
sample_2_5 <-sample_2_5[sample_2_5$GT != "0|1",]
#remove F1R2,F2R1 values that end it ",0" ensuring reads in both forward and reverse direction
sample_2_5 <- sample_2_5 %>% filter(!grepl(",0$", F1R2), !grepl(",0$", F2R1))
write.table(sample_2_5, file = "~/7965/7965-DH-0002_0005")



```

```{r}


sample_2_8 <- read.table("~/7965/7965-DH-0002_0008.vcf")
colnames(sample_2_8) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "7965-DH-0002", "7965-DH-0008")
sample_2_8 <- filter(sample_2_8, FILTER == "PASS")
#Isolate values
library(magrittr)
#GENE
GENE<-data.frame(sample_2_8[[8]]%>%strsplit("|", fixed=T)%>%sapply("[", 4))
colnames(GENE)[1] <- "GENE"
#DP
DP<-data.frame(sample_2_8[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 7))
colnames(DP)[1] <- "DP"
#F1R2
F1R2<-data.frame(sample_2_8[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 5))
colnames(F1R2)[1] <- "F1R2"
#F2R1
F2R1 <-data.frame(sample_2_8[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 6))
colnames(F2R1)[1] <- "F2R1"
#genotype
GT <-data.frame(sample_2_8[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 1))
colnames(GT)[1] <- "GT"
#combine all isolted columns to original dataframe
sample_2_8 <- cbind(sample_2_8,GENE,DP,F1R2,F2R1,GT)
sample_2_8 <-sample_2_8[,c(1:7,12,8:11,13:16)]
#DP >50
sample_2_8 <-sample_2_8[sample_2_8$DP>50,]
#remove GT = 0|1
sample_2_8 <-sample_2_8[sample_2_8$GT != "0|1",]
#remove F1R2,F2R1 values that end it ",0" ensuring reads in both forward and reverse direction
sample_2_8 <- sample_2_8 %>% filter(!grepl(",0$", F1R2), !grepl(",0$", F2R1))
write.table(sample_2_8, file = "~/7965/7965-DH-0002_0008")


```

```{r}



sample_3_6 <- read.table("~/7965/7965-DH-0003_0006.vcf")
colnames(sample_3_6) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "7965-DH-0003", "7965-DH-0006")
sample_3_6 <- filter(sample_3_6, FILTER == "PASS")
#Isolate values
library(magrittr)
#GENE
GENE<-data.frame(sample_3_6[[8]]%>%strsplit("|", fixed=T)%>%sapply("[", 4))
colnames(GENE)[1] <- "GENE"
#DP
DP<-data.frame(sample_3_6[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 7))
colnames(DP)[1] <- "DP"
#F1R2
F1R2<-data.frame(sample_3_6[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 5))
colnames(F1R2)[1] <- "F1R2"
#F2R1
F2R1 <-data.frame(sample_3_6[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 6))
colnames(F2R1)[1] <- "F2R1"
#genotype
GT <-data.frame(sample_3_6[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 1))
colnames(GT)[1] <- "GT"
#combine all isolted columns to original dataframe
sample_3_6 <- cbind(sample_3_6,GENE,DP,F1R2,F2R1,GT)
sample_3_6 <-sample_3_6[,c(1:7,12,8:11,13:16)]
#DP >50
sample_3_6 <-sample_3_6[sample_3_6$DP>50,]
#remove GT = 0|1
sample_3_6 <-sample_3_6[sample_3_6$GT != "0|1",]
#remove F1R2,F2R1 values that end it ",0" ensuring reads in both forward and reverse direction
sample_3_6 <- sample_3_6 %>% filter(!grepl(",0$", F1R2), !grepl(",0$", F2R1))
write.table(sample_3_6, file = "~/7965/7965-DH-0003_0006")


```

```{r}




sample_3_9 <- read.table("~/7965/7965-DH-0003_0009.vcf")
colnames(sample_3_9) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "7965-DH-0003", "7965-DH-0009")
sample_3_9 <- filter(sample_3_9, FILTER == "PASS")
#Isolate values
library(magrittr)
#GENE
GENE<-data.frame(sample_3_9[[8]]%>%strsplit("|", fixed=T)%>%sapply("[", 4))
colnames(GENE)[1] <- "GENE"
#DP
DP<-data.frame(sample_3_9[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 7))
colnames(DP)[1] <- "DP"
#F1R2
F1R2<-data.frame(sample_3_9[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 5))
colnames(F1R2)[1] <- "F1R2"
#F2R1
F2R1 <-data.frame(sample_3_9[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 6))
colnames(F2R1)[1] <- "F2R1"
#genotype
GT <-data.frame(sample_3_9[[11]]%>%strsplit(":", fixed=T)%>%sapply("[", 1))
colnames(GT)[1] <- "GT"
#combine all isolted columns to original dataframe
sample_3_9 <- cbind(sample_3_9,GENE,DP,F1R2,F2R1,GT)
sample_3_9 <-sample_3_9[,c(1:7,12,8:11,13:16)]
#DP >50
sample_3_9 <-sample_3_9[sample_3_9$DP>50,]
#remove GT = 0|1
sample_3_9 <-sample_3_9[sample_3_9$GT != "0|1",]
#remove F1R2,F2R1 values that end it ",0" ensuring reads in both forward and reverse direction
sample_3_9 <- sample_3_9 %>% filter(!grepl(",0$", F1R2), !grepl(",0$", F2R1))
write.table(sample_3_9, file = "~/7965/7965-DH-0003_0009")






```
