###graph of cormotif gene groups
library(limma)
library(tidyverse)
library(biomaRt)
library(RColorBrewer)
library(ggpubr)
### none work so far do eyeballing
# set1- no response_cluster24h 11190
# set3- all top2bi 1275442
# set4-late responsetop2bi only 126820
# set5-early response top2bi only 27245
drug_pal <- c("#41B333","#8B006D","#DF707E","#F1B72B", "#3386DD","#707031")
#geneexpressionsets <- cbind(sets=c('set1', 'set3','set4','set5'), ENTREZID = c(11190,127544,126820,27245))
time <- rep((rep(c("3h", "24h"), c(6,6))), 6)
time <- ordered(time, levels =c("3h", "24h"))
drug <- rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)
indv <- as.factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))
labeltop <- (interaction(substring(drug, 0, 2), indv, time))
#fills <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")

DEG_cormotif <- readRDS("data/DEG_cormotif.RDS")
motif1_NR <- DEG_cormotif$motif1_NR
motif3_TI <- DEG_cormotif$motif3_TI
motif4_LR <- DEG_cormotif$motif4_LR
motif5_ER <- DEG_cormotif$motif5_ER
set.seed(12345)# (picked 2650,7803,55588,163081)
# geneexpressionsets <- c(sample(motif1_NR,1),sample(motif3_TI,1),sample(motif4_LR,1),sample(motif5_ER,1))

geneexpressionsets=c("2650","7803","55588","27245")

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
my_chr <- c(1:22, 'M', 'X', 'Y')  ## creates a filter for each database
#attributes <- listAttributes((ensembl))
my_attributes <- c('entrezgene_id', 'ensembl_gene_id', 'hgnc_symbol')


Motif_expr<- getBM(attributes=my_attributes,filters ='entrezgene_id',
         values = geneexpressionsets, mart = ensembl)
# toplistall %>% dplyr::filter(ENTREZID==Motif_expr[1,1])

#namelist <- c("1026"="CDKN1A","23411"= "SIRT1","27113"= "BBC3","4193"= "MDM2")
cpmcounts <- read.csv("data/filtered_cpm_counts.csv", row.names = 1)
#norm_counts <- read.csv("data/norm_counts.csv", row.names = 1)
#norm_counts <- as.data.frame(log(cpm(norm_counts)))
#colnames(norm_counts) <- labeltop
#cpmcounts <- as.data.frame(cpmcounts)
colnames(cpmcounts) <- labeltop

cpmcounts %>% dplyr::filter(rownames(.)==geneexpressionsets[1]) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%
  mutate(drug=rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)) %>%
  mutate(time = rep((rep(c("3h", "24h"), c(6,6))), 6)) %>%
  mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))) %>%
  ggplot(., aes(x=drug, y=counts))+
    geom_boxplot(position="identity",aes(fill=drug))+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  guides(alpha= "none", size= "none",  indv="none")+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_pal)+
    facet_wrap("time", scales="free_y", nrow=1, ncol=2)+
  theme_bw()+
  ylab(expression(atop("No Response set",italic("GCNT1")~log[2]~"cpm ")))+
  xlab("")+
  xlab("")+
  theme(strip.background = element_rect(fill = "#C77CFF"),
        plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 12, color = "white", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))

cpmcounts %>% dplyr::filter(rownames(.)==geneexpressionsets[2]) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%
  mutate(drug=rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)) %>%
  mutate(time = rep((rep(c("3h", "24h"), c(6,6))), 6)) %>% mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))) %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot(position="identity",aes(fill=drug))+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  guides(alpha= "none", size= "none",  indv="none")+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_pal)+
  facet_wrap("time", scales="free_y", nrow=1, ncol=2)+
  theme_bw()+
  ylab(expression(atop("Time-Independent set",italic("PTP4A1")~log[2]~"cpm ")))+
  xlab("")+
  theme(strip.background = element_rect(fill = "#00BFC4"),
        plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 12, color = "white", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))

cpmcounts %>% dplyr::filter(rownames(.)==geneexpressionsets[3]) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%
  mutate(drug=rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)) %>%
  mutate(time = rep((rep(c("3h", "24h"), c(6,6))), 6)) %>%
  mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))) %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot(position="identity",aes(fill=drug))+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  guides(alpha= "none", size= "none",  indv="none")+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_pal)+
  facet_wrap("time", scales="free_y", nrow=1, ncol=2)+
  theme_bw()+
  ylab(expression(atop("Late Response set",italic("MED29")~log[2]~"cpm ")))+
  xlab("")+
  xlab("")+
  theme(strip.background = element_rect(fill = "#7CAE00"),
        plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 12, color = "white", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))

cpmcounts %>% dplyr::filter(rownames(.)==27245) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%
  mutate(drug = rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)) %>%
  mutate(time = rep((rep(c("3h", "24h"), c(6,6))), 6)) %>%
  mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))) %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot(position="identity",aes(fill=drug))+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  guides(alpha= "none", size= "none",  indv="none")+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_pal)+
  facet_wrap("time", scales="free_y", nrow=1, ncol=2)+
  theme_bw()+
  ylab(expression(atop("Early Response set",italic("AHDC1")~log[2]~"cpm ")))+
  xlab("")+
  xlab("")+
  theme(strip.background = element_rect(fill = "#F8766D"),
        plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 12, color = "white", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))

# cpm counts are above ----------------------------------------------------



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
