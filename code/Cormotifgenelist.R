###graph of cormotif gene groups

library(tidyverse)
library(biomaRt)
library(RColorBrewer)
### none work so far do eyeballing
# set1- no response_cluster24h 11190
# set3- all top2bi 1275442
# set4-late responsetop2bi only 126820
# set5-early response top2bi only 27245

#geneexpressionsets <- cbind(sets=c('set1', 'set3','set4','set5'), ENTREZID = c(11190,127544,126820,27245))
time <- rep((rep(c("3h", "24h"), c(6,6))), 6)
time <- ordered(time, levels =c("3h", "24h"))
drug <- rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)
indv <- as.factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))
labeltop <- (interaction(substring(drug, 0, 2), indv, time))


DEG_cormotif <- readRDS("data/DEG_cormotif.RDS")
motif1_NR <- DEG_cormotif$motif1_NR
motif3_TI <- DEG_cormotif$motif3_TI
motif4_LR <- DEG_cormotif$motif4_LR
motif5_ER <- DEG_cormotif$motif5_ER
set.seed(12345)
geneexpressionsets <- c(sample(motif1_NR,1),sample(motif3_TI,1),sample(motif4_LR,1),sample(motif5_ER,1))
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
my_chr <- c(1:22, 'M', 'X', 'Y')  ## creates a filter for each database
#attributes <- listAttributes((ensembl))
my_attributes <- c('entrezgene_id', 'ensembl_gene_id', 'hgnc_symbol')


Motif_expr<- getBM(attributes=my_attributes,filters ='entrezgene_id',
         values = 4356, mart = ensembl)
# toplistall %>% dplyr::filter(ENTREZID==Motif_expr[1,1])

#namelist <- c("1026"="CDKN1A","23411"= "SIRT1","27113"= "BBC3","4193"= "MDM2")

norm_counts <- read.csv("data/norm_counts.csv", row.names = 1)
norm_counts <- as.data.frame(log(cpm(norm_counts)))
colnames(norm_counts) <- labeltop

norm_counts %>% dplyr::filter(rownames(.)==geneexpressionsets[1]) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%
  mutate(drug=rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)) %>%
  mutate(time = rep((rep(c("3h", "24h"), c(6,6))), 6)) %>% mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))) %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot()+
  ylim(-1,3)+
  geom_point(aes(col=indv,size=3))+
  scale_color_brewer(palette = "Dark2",name = "Individual")+
  facet_wrap("time", scales="free_y", nrow=1, ncol=2)+
  theme_bw()+
  ylab("log cpm of ZF567 (No Response set) ")

norm_counts %>% dplyr::filter(rownames(.)==geneexpressionsets[2]) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%
  mutate(drug=rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)) %>%
  mutate(time = rep((rep(c("3h", "24h"), c(6,6))), 6)) %>% mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))) %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot()+
  geom_point(aes(col=indv,size=3))+
  scale_color_brewer(palette = "Dark2",name = "Individual")+
  facet_wrap("time",  nrow=1, ncol=2)+
  theme_bw()+
  ylab("log cpm of GCNT1 (Time independent response set) ")

norm_counts %>% dplyr::filter(rownames(.)==geneexpressionsets[3]) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%
  mutate(drug=rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)) %>%
  mutate(time = rep((rep(c("3h", "24h"), c(6,6))), 6)) %>%
  mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))) %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot()+
  geom_point(aes(col=indv,size=3))+
  scale_color_brewer(palette = "Dark2",name = "Individual")+
  facet_wrap("time", scales=, nrow=1, ncol=2)+
  theme_bw()+
  ylab("log cpm of AHDC1 (TOP2Bi-late response set) ")

norm_counts %>% dplyr::filter(rownames(.)==geneexpressionsets[4]) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%
  mutate(drug = rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)) %>%
  mutate(time = rep((rep(c("3h", "24h"), c(6,6))), 6)) %>%
  mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))) %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot()+
  geom_point(aes(col=indv,size=3))+
  scale_color_brewer(palette = "Dark2",name = "Individual")+
  facet_wrap("time",  nrow=1, ncol=2)+
  theme_bw()+
  ylab("log cpm of PTP4A1 (TOP2Bi-early response set) ")


norm_counts %>% dplyr::filter(rownames(.)==22822) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%
  mutate(drug = rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)) %>%
  mutate(time = rep((rep(c("3h", "24h"), c(6,6))), 6)) %>%
  mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))) %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot()+
  geom_point(aes(col=indv,size=3))+
  scale_color_brewer(palette = "Dark2",name = "Individual")+
  facet_wrap("time",  nrow=1, ncol=2)+
  theme_bw()+
  ylab("Motif2  PHLDA1")

norm_counts %>% dplyr::filter(rownames(.)==4356) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%
  mutate(drug = rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)) %>%
  mutate(time = rep((rep(c("3h", "24h"), c(6,6))), 6)) %>%
  mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))) %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot()+
  geom_point(aes(col=indv,size=3))+
  scale_color_brewer(palette = "Dark2",name = "Individual")+
  facet_wrap("time",  nrow=1, ncol=2)+
  theme_bw()+
  ylab("Motif MPP3")
