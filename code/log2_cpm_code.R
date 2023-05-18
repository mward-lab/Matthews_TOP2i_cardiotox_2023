library(tidyverse)
library(limma)

library(biomaRt)
library(RColorBrewer)
library(ggpubr)

### regex to trim things after a name     mutate(name= sub("\\.\\d+$", "", name)) %>%
DEG_cormotif <- readRDS("data/DEG_cormotif.RDS")
list2env(DEG_cormotif,envir=.GlobalEnv)


toplistall <- readRDS("data/toplistall.RDS")

fil_list_gene <- c("RARG", "TOP2B")
fil_list_Entrez <- c("5916", "7155")
toplistall %>%
  #filter(time =="3_hours") %>%
  filter(SYMBOL=="RARG"| SYMBOL=="TOP2B") %>%
  ggplot(., aes(x=id, y=logFC, color=SYMBOL))+
  geom_point()+
  ylim(-2,2)+
  facet_wrap(~time)

cpmcounts <- read.csv("data/filtered_cpm_counts.csv", row.names = 1)

efit2<- readRDS(file ="data/efit2results.RDS")
  time_lb <- rep((rep(c("3h", "24h"), c(6,6))), 6)
  time_lb <- ordered(time_lb, levels =c("3h", "24h"))
  drug_lb <- rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Vehicle"),12)
  indv_lb <- as.factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))
  labeltop_noInd <- (interaction(substring(drug_lb, 0, 2), time_lb))
cpmcounts %>%
  set_names(labeltop) %>%
  rownames_to_column() %>%
  filter(rowname=="5916"| rowname=="7155") %>%
  pivot_longer(!rowname, names_to = "treatment",values_to = "counts") %>%
  separate(treatment,c("drug","indv","time")) %>%
  mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(indv, levels = c(1,2,3,4,5,6))) %>%
  mutate(drug =case_match(drug, "Da"~"Daunorubicin",
                          "Do"~"Doxorubicin",
                          "Ep"~"Epirubicin",
                          "Mi"~"Mitoxantrone",
                          "Tr"~"Trastuzumab",
                          "Ve"~"Vehicle", .default = drug)) %>%
  filter(time =="24h") %>%
  ggplot(., aes(x=drug,  y=counts))+
  geom_boxplot()+
  facet_wrap(~rowname)



# CDKN1a and Top2B log2CPM ------------------------------------------------


drug_pal_vehend <- c("#8B006D","#DF707E","#F1B72B", "#3386DD","#707031","#41B333")

time <- rep((rep(c("3h", "24h"), c(6,6))), 6)
time <- ordered(time, levels =c("3h", "24h"))
drug <- rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Vehicle"),12)
indv <- as.factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))
labeltop <- (interaction(substring(drug, 0, 2), indv, time))


ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
my_chr <- c(1:22, 'M', 'X', 'Y')  ## creates a filter for each database

my_attributes <- c('entrezgene_id', 'ensembl_gene_id', 'hgnc_symbol')


cpmcounts <- read.csv("data/filtered_cpm_counts.csv", row.names = 1)

### TOP2B
cpmcounts %>% dplyr::filter(rownames(.)==7155) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%

  separate(treatment,c("drug","indv","time")) %>%

  mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(indv, levels = c(1,2,3,4,5,6))) %>%
  mutate(drug =case_match(drug, "Da"~"Daunorubicin",
                          "Do"~"Doxorubicin",
                          "Ep"~"Epirubicin",
                          "Mi"~"Mitoxantrone",
                          "Tr"~"Trastuzumab",
                          "Ve"~"Vehicle", .default = drug))  %>%

  filter(time=="24h") %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot(position="identity",aes(fill=drug))+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  guides(alpha= "none", size= "none")+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_pal_vehend)+
  geom_signif(comparisons =list(c("Daunorubicin","Vehicle"),
                                c("Doxorubicin","Vehicle"),
                                c("Epirubicin","Vehicle"),
                                c("Mitoxantrone","Vehicle"),
                                c("Trastuzumab","Vehicle")),
              test= "t.test",
              map_signif_level=FALSE,
              textsize =6,
              tip_length = .1,
              vjust = 0.2,
              step_increase = 0.1)+
  facet_wrap("time", scales="free_y", nrow=1, ncol=2)+
  theme_bw()+
  ylab(expression(atop(" ",italic("Top2"*beta)~log[2]~"cpm ")))+
  xlab("")+
  xlab("")+
  theme(#strip.background = element_rect(fill = "#C77CFF"),
        plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 12, color = "white", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))

### CDKN1a
cpmcounts %>% dplyr::filter(rownames(.)==1026) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%

  separate(treatment,c("drug","indv","time")) %>%

  mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(indv, levels = c(1,2,3,4,5,6))) %>%
  mutate(drug =case_match(drug, "Da"~"Daunorubicin",
                          "Do"~"Doxorubicin",
                          "Ep"~"Epirubicin",
                          "Mi"~"Mitoxantrone",
                          "Tr"~"Trastuzumab",
                          "Ve"~"Vehicle", .default = drug))  %>%

  filter(time=="24h") %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot(position="identity",aes(fill=drug))+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  guides(alpha= "none", size= "none")+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_pal_vehend)+
  geom_signif(comparisons =list(c("Daunorubicin","Vehicle"),
                                c("Doxorubicin","Vehicle"),
                                c("Epirubicin","Vehicle"),
                                c("Mitoxantrone","Vehicle"),
                                c("Trastuzumab","Vehicle")),
              test= "t.test",
              map_signif_level=FALSE,
              textsize =6,
              tip_length = .1,
              vjust = 0.2,
              step_increase = 0.1)+
  facet_wrap("time", scales="free_y", nrow=1, ncol=2)+
  theme_bw()+
  ylab(expression(atop(" ",italic("CDKN1"*alpha)~log[2]~"cpm ")))+
  xlab("")+
  xlab("")+
  theme(strip.background = element_rect(fill = "#7CAE00"),
        plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 12, color = "white", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))


### Top2a

cpmcounts %>% dplyr::filter(rownames(.)==7153) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%

  separate(treatment,c("drug","indv","time")) %>%

  mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(indv, levels = c(1,2,3,4,5,6))) %>%
  mutate(drug =case_match(drug, "Da"~"Daunorubicin",
                          "Do"~"Doxorubicin",
                          "Ep"~"Epirubicin",
                          "Mi"~"Mitoxantrone",
                          "Tr"~"Trastuzumab",
                          "Ve"~"Vehicle", .default = drug))  %>%

  filter(time=="24h") %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot(position="identity",aes(fill=drug))+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  guides(alpha= "none", size= "none")+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_pal_vehend)+
  geom_signif(comparisons =list(c("Daunorubicin","Vehicle"),
                                c("Doxorubicin","Vehicle"),
                                c("Epirubicin","Vehicle"),
                                c("Mitoxantrone","Vehicle"),
                                c("Trastuzumab","Vehicle")),
              test= "t.test",
              map_signif_level=FALSE,
              textsize =6,
              tip_length = .1,
              vjust = 0.2,
              step_increase = 0.1)+
  facet_wrap("time", scales="free_y", nrow=1, ncol=2)+
  theme_bw()+
  ylab(expression(atop(" ",italic("Top2"*alpha)~log[2]~"cpm ")))+
  xlab("")+
  xlab("")+
  theme(strip.background = element_rect(fill = "white"),
    plot.title = element_text(size=18,hjust = 0.5),
    axis.title = element_text(size = 15, color = "black"),
    axis.ticks = element_line(linewidth = 1.5),
    axis.line = element_line(linewidth = 1.5),
    axis.text.x = element_text(size = 12, color = "white", angle = 0),
    strip.text.x = element_text(size = 15, color = "black", face = "bold"))
### ATM

cpmcounts %>% dplyr::filter(rownames(.)==472) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%

  separate(treatment,c("drug","indv","time")) %>%

  mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(indv, levels = c(1,2,3,4,5,6))) %>%
  mutate(drug =case_match(drug, "Da"~"Daunorubicin",
                          "Do"~"Doxorubicin",
                          "Ep"~"Epirubicin",
                          "Mi"~"Mitoxantrone",
                          "Tr"~"Trastuzumab",
                          "Ve"~"Vehicle", .default = drug))  %>%

  filter(time=="24h") %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot(position="identity",aes(fill=drug))+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  guides(alpha= "none", size= "none")+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_pal_vehend)+
  geom_signif(comparisons =list(c("Daunorubicin","Vehicle"),
                                c("Doxorubicin","Vehicle"),
                                c("Epirubicin","Vehicle"),
                                c("Mitoxantrone","Vehicle"),
                                c("Trastuzumab","Vehicle")),
              test= "t.test",
              map_signif_level=FALSE,
              textsize =6,
              tip_length = .1,
              vjust = 0.2,
              step_increase = 0.1)+
  facet_wrap("time", scales="free_y", nrow=1, ncol=2)+
  theme_bw()+
  ylab(expression(atop(" ",italic("ATM")~log[2]~"cpm ")))+
  xlab("")+
  xlab("")+
  theme(#strip.background = element_rect(fill = "#C77CFF"),
    plot.title = element_text(size=18,hjust = 0.5),
    axis.title = element_text(size = 15, color = "black"),
    axis.ticks = element_line(linewidth = 1.5),
    axis.line = element_line(linewidth = 1.5),
    axis.text.x = element_text(size = 12, color = "white", angle = 0),
    strip.text.x = element_text(size = 15, color = "black", face = "bold"))

###  ATR

cpmcounts %>% dplyr::filter(rownames(.)==545) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%

  separate(treatment,c("drug","indv","time")) %>%

  mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(indv, levels = c(1,2,3,4,5,6))) %>%
  mutate(drug =case_match(drug, "Da"~"Daunorubicin",
                          "Do"~"Doxorubicin",
                          "Ep"~"Epirubicin",
                          "Mi"~"Mitoxantrone",
                          "Tr"~"Trastuzumab",
                          "Ve"~"Vehicle", .default = drug))  %>%

  filter(time=="24h") %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot(position="identity",aes(fill=drug))+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  guides(alpha= "none", size= "none")+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_pal_vehend)+
  geom_signif(comparisons =list(c("Daunorubicin","Vehicle"),
                                c("Doxorubicin","Vehicle"),
                                c("Epirubicin","Vehicle"),
                                c("Mitoxantrone","Vehicle"),
                                c("Trastuzumab","Vehicle")),
              test= "t.test",
              map_signif_level=FALSE,
              textsize =6,
              tip_length = .1,
              vjust = 0.2,
              step_increase = 0.1)+
  facet_wrap("time", scales="free_y", nrow=1, ncol=2)+
  theme_bw()+
  ylab(expression(atop(" ",italic("ATR")~log[2]~"cpm ")))+
  xlab("")+
  xlab("")+
  theme(#strip.background = element_rect(fill = "#C77CFF"),
    plot.title = element_text(size=18,hjust = 0.5),
    axis.title = element_text(size = 15, color = "black"),
    axis.ticks = element_line(linewidth = 1.5),
    axis.line = element_line(linewidth = 1.5),
    axis.text.x = element_text(size = 12, color = "white", angle = 0),
    strip.text.x = element_text(size = 15, color = "black", face = "bold"))

##Rictor


cpmcounts %>% dplyr::filter(rownames(.)==545) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%

  separate(treatment,c("drug","indv","time")) %>%

  mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(indv, levels = c(1,2,3,4,5,6))) %>%
  mutate(drug =case_match(drug, "Da"~"Daunorubicin",
                          "Do"~"Doxorubicin",
                          "Ep"~"Epirubicin",
                          "Mi"~"Mitoxantrone",
                          "Tr"~"Trastuzumab",
                          "Ve"~"Vehicle", .default = drug))  %>%

  filter(time=="24h") %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot(position="identity",aes(fill=drug))+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  guides(alpha= "none", size= "none")+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_pal_vehend)+
  geom_signif(comparisons =list(c("Daunorubicin","Vehicle"),
                                c("Doxorubicin","Vehicle"),
                                c("Epirubicin","Vehicle"),
                                c("Mitoxantrone","Vehicle"),
                                c("Trastuzumab","Vehicle")),
              test= "t.test",
              map_signif_level=FALSE,
              textsize =6,
              tip_length = .1,
              vjust = 0.2,
              step_increase = 0.1)+
  facet_wrap("time", scales="free_y", nrow=1, ncol=2)+
  theme_bw()+
  ylab(expression(atop(" ",italic("RICTOR")~log[2]~"cpm ")))+
  xlab("")+
  xlab("")+
  theme(#strip.background = element_rect(fill = "#C77CFF"),
    plot.title = element_text(size=18,hjust = 0.5),
    axis.title = element_text(size = 15, color = "black"),
    axis.ticks = element_line(linewidth = 1.5),
    axis.line = element_line(linewidth = 1.5),
    axis.text.x = element_text(size = 12, color = "white", angle = 0),
    strip.text.x = element_text(size = 15, color = "black", face = "bold"))

### mTOR


cpmcounts %>% dplyr::filter(rownames(.)==2475) %>%
  pivot_longer(everything(), names_to = "treatment",values_to = "counts") %>%

  separate(treatment,c("drug","indv","time")) %>%

  mutate(time=factor(time, levels =c("3h", "24h"))) %>%
  mutate(indv=factor(indv, levels = c(1,2,3,4,5,6))) %>%
  mutate(drug =case_match(drug, "Da"~"Daunorubicin",
                          "Do"~"Doxorubicin",
                          "Ep"~"Epirubicin",
                          "Mi"~"Mitoxantrone",
                          "Tr"~"Trastuzumab",
                          "Ve"~"Vehicle", .default = drug))  %>%

  filter(time=="24h") %>%
  ggplot(., aes(x=drug, y=counts))+
  geom_boxplot(position="identity",aes(fill=drug))+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  guides(alpha= "none", size= "none")+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_pal_vehend)+
  geom_signif(comparisons =list(c("Daunorubicin","Vehicle"),
                                c("Doxorubicin","Vehicle"),
                                c("Epirubicin","Vehicle"),
                                c("Mitoxantrone","Vehicle"),
                                c("Trastuzumab","Vehicle")),
              test= "t.test",
              map_signif_level=FALSE,
              textsize =6,
              tip_length = .1,
              vjust = 0.2,
              step_increase = 0.1)+
  facet_wrap("time", scales="free_y", nrow=1, ncol=2)+
  theme_bw()+
  ylab(expression(atop(" ",italic("mTOR")~log[2]~"cpm ")))+
  xlab("")+
  xlab("")+
  theme(#strip.background = element_rect(fill = "#C77CFF"),
    plot.title = element_text(size=18,hjust = 0.5),
    axis.title = element_text(size = 15, color = "black"),
    axis.ticks = element_line(linewidth = 1.5),
    axis.line = element_line(linewidth = 1.5),
    axis.text.x = element_text(size = 12, color = "white", angle = 0),
    strip.text.x = element_text(size = 15, color = "black", face = "bold"))


