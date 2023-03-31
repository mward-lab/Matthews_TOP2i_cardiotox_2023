###analysis of Eqtl from Knowles paper



# libraries ---------------------------------------------------------------


library(limma)
library(tidyverse)
library(ggsignif)
library(biomaRt)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(scales)










# loading files/counts,etc ------------------------------------------------


drug_palc <- c("#8B006D","#DF707E","#F1B72B", "#3386DD","#707031","#41B333")
time <- rep((rep(c("3h", "24h"), c(6,6))), 6)
time <- ordered(time, levels =c("3h", "24h"))
drug <- rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)
indv <- as.factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))
labeltop <- (interaction(substring(drug, 0, 2), indv, time))

level_order2 <- c('75','87','77','79','78','71')


efit2<- readRDS(file ="data/efit2results.RDS")
V.DA.top= topTable(efit2, coef=1, adjust="BH", number=Inf, sort.by="p")
V.DX.top= topTable(efit2, coef=2, adjust="BH", number=Inf, sort.by="p")
V.EP.top= topTable(efit2, coef=3, adjust="BH", number=Inf, sort.by="p")
V.MT.top= topTable(efit2, coef=4, adjust="BH", number=Inf, sort.by="p")
V.TR.top= topTable(efit2, coef=5, adjust="BH", number=Inf, sort.by="p")

toplist3hours <-list(V.DA.top[,c(1:3,4,6:7)],
                     V.DX.top[,c(1:3,4,6:7)],
                     V.EP.top[,c(1:3,4,6:7)],
                     V.MT.top[,c(1:3,4,6:7)],
                     V.TR.top[,c(1:3,4,6:7)])
names(toplist3hours) <- c("Daunorubicin","Doxorubicin", "Epirubicin","Mitoxantrone", "Trastuzumab")

toplist3hours <- map_df(toplist3hours, ~as.data.frame(.x), .id="id")



V.DA24.top= topTable(efit2, coef=6, adjust="BH", number=Inf, sort.by="p")
V.DX24.top= topTable(efit2, coef=7, adjust="BH", number=Inf, sort.by="p")
V.EP24.top= topTable(efit2, coef=8, adjust="BH", number=Inf, sort.by="p")
V.MT24.top= topTable(efit2, coef=9, adjust="BH", number=Inf, sort.by="p")
V.TR24.top= topTable(efit2, coef=10, adjust="BH", number=Inf, sort.by="p")


toplist24hours <-list(V.DA24.top[,c(1:3,4,6:7)],
                      V.DX24.top[,c(1:3,4,6:7)],
                      V.EP24.top[,c(1:3,4,6:7)],
                      V.MT24.top[,c(1:3,4,6:7)],
                      V.TR24.top[,c(1:3,4,6:7)])
names(toplist24hours) <- c("Daunorubicin","Doxorubicin", "Epirubicin","Mitoxantrone", "Trastuzumab")

toplist24hours <- map_df(toplist24hours, ~as.data.frame(.x), .id="id")

toplistall <- list(toplist24hours,toplist3hours)
names(toplistall) <- c("24_hours", "3_hours")
toplistall <-map_df(toplistall, ~as.data.frame(.x), .id="time")


#-Knowles data:
#fig4  contains all sig marginal effect eqtl
#fig5 contains e1gl using total expression only genes
#fig6 contains eqtl mapped using total and allele specificing expression

knowfig5<- Knowles_2018.elife.33480.supp5.v2 <- read.delim("~/Ward Lab/Cardiotoxicity/Manuscript/Knowles_2018-elife-33480-supp5-v2/Knowles_2018-elife-33480-supp5-v2")

knowfig6 <- Knowles_2018.elife.33480.supp6.v2 <- read.delim("~/Ward Lab/Cardiotoxicity/Manuscript/Knowles_2018-elife-33480-supp6-v2/Knowles_2018-elife-33480-supp6-v2")

knowfig4 <- Knowles_2018.elife.33480.supp4.v2 <- read.delim("~/Ward Lab/Cardiotoxicity/Manuscript/Knowles_2018-elife-33480-supp4-v2/Knowles_2018-elife-33480-supp4-v2")


# add in ensemble database ------------------------------------------------

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
my_chr <- c(1:22, 'M', 'X', 'Y')  ## creates a filter for each
attributes <- listAttributes((ensembl))
my_attributes <- c('entrezgene_id', 'ensembl_gene_id', 'hgnc_symbol')


# covert ensg to entrezid -------------------------------------------------

knowles5 <- getBM(attributes=my_attributes,filters ='ensembl_gene_id',
                  values =knowfig5, mart = ensembl)
#377
knowles6 <- getBM(attributes=my_attributes,filters ='ensembl_gene_id',
                  values =knowfig6, mart = ensembl)
#450
knowles4 <- getBM(attributes=my_attributes,filters ='ensembl_gene_id',
                  values =knowfig4, mart = ensembl)
#524

knowles4 <-unique(knowles4$entrezgene_id)
knowles5 <-unique(knowles5$entrezgene_id)
knowles6 <-unique(knowles6$entrezgene_id)

all_QTLSknow<- (union(union(knowles5,knowles6),knowles4))




#
# toplistall %>% subset(ENTREZID %in% knowles4)



# knowles4 ----------------------------------------------------------------



toplistall %>% filter(adj.P.Val<0.1) %>%
  filter(ENTREZID %in% knowles4) %>%
  mutate(id = as.factor(id)) %>%
  mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%
  ggplot(., aes(x=id, y =logFC))+
  geom_boxplot(aes(fill=id))+
  fill_palette(palette =drug_palNoVeh)+
  guides(fill=guide_legend(title = "Treatment"))+
  #facet_wrap(~time,labeller = (time = facettimelabel) )+
  theme_bw()+
  xlab("")+
  ylab("Log Fold Change")+
  theme_bw()+
  facet_wrap(~time)+
  ggtitle("all expressed eQTL (supp 4)")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        # axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        strip.background = element_rect(fill = "transparent"),
        axis.text = element_text(size = 8, color = "black", angle = 0),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  geom_signif(comparisons = list(c("Daunorubicin", "Doxorubicin"),
                                 c("Epirubicin", "Doxorubicin"),
                                 c("Mitoxantrone", "Doxorubicin")),
                              test = "t.test",
              map_signif_level = FALSE,step_increase = 0.1,
              textsize = 4)

# knowles5 ----------------------------------------------------------------

toplistall %>%
  filter(adj.P.Val<0.1) %>%
  filter(ENTREZID %in% knowles5) %>%
  mutate(id = as.factor(id)) %>%
  mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%
  ggplot(., aes(x=id, y =logFC))+
  geom_boxplot(aes(fill=id))+
  geom_point()+
  fill_palette(palette =drug_palNoVeh)+
  guides(fill=guide_legend(title = "Treatment"))+
  #facet_wrap(~time,labeller = (time = facettimelabel) )+
  theme_bw()+
  xlab("")+
  ylab("Log Fold Change")+
  theme_bw()+
  facet_wrap(~time)+
  ggtitle("Response eQTL (supp 5)")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        # axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        strip.background = element_rect(fill = "transparent"),
        axis.text = element_text(size = 8, color = "black", angle = 0),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+


  geom_signif(comparisons = list(c("Daunorubicin", "Doxorubicin"),
                                 c("Epirubicin", "Doxorubicin"),
                                 c("Mitoxantrone", "Doxorubicin")),

              test = "t.test",
              map_signif_level = FALSE,step_increase = 0.1,
              textsize = 4)

# knowles 6 ---------------------------------------------------------------

toplistall %>%filter(adj.P.Val<0.1) %>%
  filter(ENTREZID %in% knowles6) %>%
  mutate(id = as.factor(id)) %>%
  mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%
  ggplot(., aes(x=id, y =logFC))+
  geom_boxplot(aes(fill=id))+
  fill_palette(palette =drug_palNoVeh)+
  guides(fill=guide_legend(title = "Treatment"))+
  #facet_wrap(~time,labeller = (time = facettimelabel) )+
  theme_bw()+
  xlab("")+
  ylab("Log Fold Change")+
  theme_bw()+
  facet_wrap(~time)+
  ggtitle("splice eQTL (supp 6)")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        # axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        strip.background = element_rect(fill = "transparent"),
        axis.text = element_text(size = 8, color = "black", angle = 0),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+


  geom_signif(comparisons = list(c("Daunorubicin", "Doxorubicin"),

                                 c("Epirubicin", "Doxorubicin"),
                                 c("Mitoxantrone", "Doxorubicin")),

              test = "t.test",
              map_signif_level = FALSE,step_increase = 0.1,
              textsize = 4)

# union knowles all -----------------------------------
filter(toplistall$ENTREZID %in% union)

all_QTLSknow<- (union(union(knowles5,knowles6),knowles4))
non_eQTL <- setdiff(toplistall$ENTREZID, all_QTLSknow)


toplistall %>%
  dplyr::filter(as.numeric(adj.P.Val)<0.1) %>%
  filter(ENTREZID %in% non_eQTL) %>%
  mutate(id = as.factor(id)) %>%
  mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%
  slice_max(ENTREZID>100000) %>%
  ggplot(., aes(x=id, y =logFC))+
  geom_boxplot(aes(fill=id))+
  fill_palette(palette =drug_palNoVeh)+
  guides(fill=guide_legend(title = "Treatment"))+
  #facet_wrap(~time,labeller = (time = facettimelabel) )+
  theme_bw()+
  xlab("")+
  ylab("Log Fold Change")+
  theme_bw()+
  facet_wrap(ENTREZID~time)+
  ggtitle("non-eQTL sigDEGs")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        # axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        strip.background = element_rect(fill = "transparent"),
        axis.text = element_text(size = 8, color = "black", angle = 0),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))#+


  geom_signif(comparisons = list(c("Daunorubicin", "Doxorubicin"),

                                 c("Epirubicin", "Doxorubicin"),
                                 c("Mitoxantrone", "Doxorubicin")),

              test = "t.test",
              map_signif_level = FALSE,step_increase = 0.1,
              textsize = 4)


sigtoplistall = toplistall[toplistall$adj.P.Val < .1 , ]
  sigtoplistall %>%
    filter(ENTREZID %in% non_eQTL)
  #%>% group_by(time,symbol) %>%

   count(ENTREZID, sort=TRUE)

 length(unique(sigtoplistall$ENTREZID))
#9047
