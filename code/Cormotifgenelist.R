###graph of cormotif gene groups

library(tidyverse)
library(biomaRt)
library(RColorBrewer)
### none work so far do eyeballing
# set1- no response_cluster24h 11190
# set3- all top2bi 1275442
# set4-late responsetop2bi only 126820
# set5-early response top2bi only 27245

geneexpressionsets <- cbind(sets=c('set1', 'set3','set4','set5'), ENTREZID = c(11190,127544,126820,27245))
time <- rep((rep(c("3h", "24h"), c(6,6))), 6)
time <- ordered(time, levels =c("3h", "24h"))
drug <- rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)
indv <- as.factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))
labeltop <- (interaction(substring(drug, 0, 2), indv, time))




ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
my_chr <- c(1:22, 'M', 'X', 'Y')  ## creates a filter for each database
#attributes <- listAttributes((ensembl))
my_attributes <- c('entrezgene_id', 'ensembl_gene_id', 'hgnc_symbol')


Motif_expr<- getBM(attributes=my_attributes,filters ='entrezgene_id',
         values = geneexpressionsets, mart = ensembl)
# toplistall %>% dplyr::filter(ENTREZID==Motif_expr[1,1])

#namelist <- c("1026"="CDKN1A","23411"= "SIRT1","27113"= "BBC3","4193"= "MDM2")

norm_counts <- read.csv("data/norm_counts.csv", row.names = 1)

colnames(norm_counts) <- labeltop

norm_counts %>% dplyr::filter(rownames(.)==Motif_expr$entrezgene_id[1]) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%
  mutate(drug=rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)) %>%
  mutate(time = rep((rep(c("3h", "24h"), c(6,6))), 6)) %>% mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))) %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot()+
  geom_point(aes(col=indv,size=3))+
  scale_color_brewer(palette = "Dark2",name = "Individual")+
  facet_wrap("time", scales="free_y", nrow=1, ncol=2)+
  theme_bw()+
  ylab("log cpm of CENP250 (No Response set) ")

norm_counts %>% dplyr::filter(rownames(.)==Motif_expr$entrezgene_id[2]) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%
  mutate(drug=rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)) %>%
  mutate(time = rep((rep(c("3h", "24h"), c(6,6))), 6)) %>% mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))) %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot()+
  geom_point(aes(col=indv,size=3))+
  scale_color_brewer(palette = "Dark2",name = "Individual")+
  facet_wrap("time", scales="free_y", nrow=1, ncol=2)+
  theme_bw()+
  ylab("log cpm of DNAI3 (Late TOP2Bi response set) ")

norm_counts %>% dplyr::filter(rownames(.)==Motif_expr$entrezgene_id[4]) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%
  mutate(drug=rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)) %>%
  mutate(time = rep((rep(c("3h", "24h"), c(6,6))), 6)) %>% mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))) %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot()+
  geom_point(aes(col=indv,size=3))+
  scale_color_brewer(palette = "Dark2",name = "Individual")+
  facet_wrap("time", scales="free_y", nrow=1, ncol=2)+
  theme_bw()+
  ylab("log cpm of AHDC1 (TOP2Bi-early response set) ")




group_by(Samples) %>%
  ggplot(.,aes(x= Samples, y=Counts))+
  geom_boxplot(aes(col=Samples))
