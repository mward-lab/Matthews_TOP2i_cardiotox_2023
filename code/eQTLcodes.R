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
library(sjmisc)
# loading files/counts,etc ------------------------------------------------

drug_palNoVeh <- c("#8B006D","#DF707E","#F1B72B", "#3386DD","#707031")
drug_palc <- c("#8B006D","#DF707E","#F1B72B", "#3386DD","#707031","#41B333")
time <- rep((rep(c("3h", "24h"), c(6,6))), 6)
time <- ordered(time, levels =c("3h", "24h"))
drug <- rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Vehicle"),12)
indv <- as.factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))
labeltop <- (interaction(substring(drug, 0, 2), indv, time))

level_order2 <- c('75','87','77','79','78','71')


efit2<- readRDS(file ="data/efit2.RDS")
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
#write.csv(toplistall, "output/toplistall.csv")


#-Knowles data:
#fig4  contains all sig marginal effect eqtl
#fig5 contains e1gl using total expression only genes
#fig6 contains eqtl mapped using total and allele specific expression not using anymore

knowfig5<- Knowles_2018.elife.33480.supp5.v2 <- read.delim("~/Ward Lab/Cardiotoxicity/Manuscript/Knowles_2018-elife-33480-supp5-v2/Knowles_2018-elife-33480-supp5-v2")

#knowfig6 <- Knowles_2018.elife.33480.supp6.v2 <- read.delim("~/Ward Lab/Cardiotoxicity/Manuscript/Knowles_2018-elife-33480-supp6-v2/Knowles_2018-elife-33480-supp6-v2")

knowfig4 <- Knowles_2018.elife.33480.supp4.v2 <- read.delim("~/Ward Lab/Cardiotoxicity/Manuscript/Knowles_2018-elife-33480-supp4-v2/Knowles_2018-elife-33480-supp4-v2")

file.names <- list.files(path = "data/", pattern = "sig*", ignore.case = TRUE,full.names = TRUE)
file.names <- file.names[c(2:10)]


filenameonly <- read_csv("data/filenameonly.txt")


for (k in 1:length(file.names)){

  assign(paste0(filenameonly$x[k]) , read.csv(file.names[k]))
}

colnames(sigVDA24)<- c("ENTREZID","SYMBOL","logFC","AveExpr","t","P.Value","adj.P.Val","B")
colnames(sigVDX24)<- c("ENTREZID","SYMBOL","logFC","AveExpr","t","P.Value","adj.P.Val","B")
colnames(sigVEP24)<- c("ENTREZID","SYMBOL","logFC","AveExpr","t","P.Value","adj.P.Val","B")
colnames(sigVMT24)<- c("ENTREZID","SYMBOL","logFC","AveExpr","t","P.Value","adj.P.Val","B")
colnames(sigVTR24)<- c("ENTREZID","SYMBOL","logFC","AveExpr","t","P.Value","adj.P.Val","B")
colnames(sigVDA3)<- c("ENTREZID","SYMBOL","logFC","AveExpr","t","P.Value","adj.P.Val","B")
colnames(sigVDX3)<- c("ENTREZID","SYMBOL","logFC","AveExpr","t","P.Value","adj.P.Val","B")
colnames(sigVEP3)<- c("ENTREZID","SYMBOL","logFC","AveExpr","t","P.Value","adj.P.Val","B")
colnames(sigVMT3)<- c("ENTREZID","SYMBOL","logFC","AveExpr","t","P.Value","adj.P.Val","B")
colnames(sigVTR3)<- c("ENTREZID","SYMBOL","logFC","AveExpr","t","P.Value","adj.P.Val","B")
backGL <- read.csv("data/backGL.txt")


# add in ensemble database ------------------------------------------------

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
my_chr <- c(1:22, 'M', 'X', 'Y')  ## creates a filter for each
attributes <- listAttributes((ensembl))
my_attributes <- c('entrezgene_id', 'ensembl_gene_id', 'hgnc_symbol')


# covert ensg to entrezid -------------------------------------------------

knowles5 <- getBM(attributes=my_attributes,filters ='ensembl_gene_id',
                  values =knowfig5, mart = ensembl)
#377
# knowles6 <- getBM(attributes=my_attributes,filters ='ensembl_gene_id',
#                   values =knowfig6, mart = ensembl)
# #450
knowles4 <- getBM(attributes=my_attributes,filters ='ensembl_gene_id',
                  values =knowfig4, mart = ensembl)
#524

knowles4 <-unique(knowles4$entrezgene_id)
knowles5 <-unique(knowles5$entrezgene_id)
# knowles6 <-unique(knowles6$entrezgene_id)

all_QTLSknow<- union(knowles5$entrezgene_id,knowles4$entrezgene_id)


#total 794


# sig_DEGname.noqtl <- setdiff(all_sigDEGname,all_QTLSknow)#8260
# intersect(all_sigDEG, all_QTLSknow)
# intersect(all_sigDEG,knowles4)#259
# listK4 <- intersect(all_sigDEG,knowles4)
# knowles4 ----------------------------------------------------------------
length(intersect(knowles4,sigVDA24$ENTREZID ))
#239
length(intersect(knowles4,sigVDA3$ENTREZID ))
#8
length(intersect(knowles4,sigVDA24$ENTREZID ))

length(intersect(knowles4,sigVDA24$ENTREZID ))

 length(intersect(knowles4,sigVDA24$ENTREZID ))


length(intersect(knowles4,test))
#503/521
test <- backGL$ENTREZID

allset <- (toplistall[toplistall$ENTREZID%in% knowles4,3])
unique(allset)

toplistall %>%
  mutate(id = as.factor(id)) %>%
  mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%
  filter(time== "24_hours") %>%

  group_by(time, id) %>%

  mutate(K4 = if_else(ENTREZID %in% knowles4,1,0))%>%
  mutate(K5 = if_else(ENTREZID %in% knowles5,1,0))%>%
  filter(adj.P.Val<0.1) %>%
 summarize(n=n())



  # count(sigcount) %>%
  # pivot_wider(id_cols = c(time,id), names_from=sigcount, values_from=n) %>%
  # mutate(prop = sig/(sig+notsig)*100) %>%
  # mutate(prop=if_else(is.na(prop),0,prop)) %>%
  ggplot(., aes(x=id ,y=K4, group=id))+
  geom_col()#+
#geom_stat(sum(K4=1))
  geom_text(aes(label = sprintf("%.2f",prop)), position=position_dodge(0.9),vjust=-.2 )+
  scale_fill_manual(values =drug_palNoVeh)
  # filter(adj.P.Val<0.1) %>%

  ggplot(., aes(x=id, y =abs(logFC)))+
  geom_boxplot(aes(fill=id))+
  fill_palette(palette =drug_palNoVeh)+
  guides(fill=guide_legend(title = "Treatment"))+
  #facet_wrap(~time,labeller = (time = facettimelabel) )+
  theme_bw()+
  xlab("")+
  ylab("abs |Log Fold Change|")+
  theme_bw()+
  facet_wrap(~time)+
  ggtitle("All QTLs from all expressed genes (n=507)")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        # axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        strip.background = element_rect(fill = "transparent"),
        axis.text = element_text(size = 8, color = "black", angle = 0),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  geom_signif(comparisons = list(c("Daunorubicin", "Doxorubicin"),
                                 c("Epirubicin", "Doxorubicin"),
                                 c("Mitoxantrone", "Doxorubicin"),
                                   c("Trastuzumab","Doxorubicin"),
                                 c("Mitoxantrone","Trastuzumab")),
                              test = "t.test",
              map_signif_level = FALSE,step_increase = 0.1,
              textsize = 4)

# knowles5 ----------------------------------------------------------------
intersect(test,knowles5)#374/377
listK5 <- intersect(all_sigDEG,knowles5)
toplistall %>%
  # filter(adj.P.Val<0.1) %>%
  filter(ENTREZID %in% listK5) %>%
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


allset <- (toplistall[toplistall$ENTREZID%in% knowles5,3])
unique(allset)


toplistall %>%
  filter(ENTREZID %in% knowles5) %>%
  mutate(id = as.factor(id)) %>%
  mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%
  ggplot(., aes(x=id, y =abs(logFC)))+
  geom_boxplot(aes(fill=id))+
  fill_palette(palette =drug_palNoVeh)+
  guides(fill=guide_legend(title = "Treatment"))+
  #facet_wrap(~time,labeller = (time = facettimelabel) )+
  theme_bw()+
  xlab("")+
  ylab("abs |Log Fold Change|")+
  theme_bw()+
  facet_wrap(~time)+
  ggtitle("All reQTLs from all expressed genes (n=375)")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        # axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        strip.background = element_rect(fill = "transparent"),
        axis.text = element_text(size = 8, color = "black", angle = 0),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  geom_signif(comparisons = list(c("Daunorubicin", "Doxorubicin"),
                                 c("Epirubicin", "Doxorubicin"),
                                 c("Mitoxantrone", "Doxorubicin"),
                                 c("Trastuzumab","Doxorubicin")),
              test = "t.test",
              map_signif_level = FALSE,step_increase = 0.1,
              textsize = 4)

# knowles 6 ---------------------------------------------------------------
intersect(all_sigDEG,knowles6)#234
listK6 <- intersect(all_sigDEG,knowles6)
toplistall %>%
  # filter(adj.P.Val<0.1) %>%
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

all_QTLSknow<- (union(knowles5,knowles4)) #794
sig3all<- union(sigVMT3$ENTREZID,union(sigVEP3$ENTREZID,union(sigVDA3$ENTREZID,sigVDX3$ENTREZID)))
sig24all<- union(sigVMT24$ENTREZID,union(sigVEP24$ENTREZID,union(sigVDA24$ENTREZID,sigVDX24$ENTREZID)))
length(unique(sig24all))
sig_3_24_all <- union(sig3all,sig24all)  #9491
NR_respdeg <- setdiff(backGL$ENTREZID, sig_3_24_all)  #4593


non_eQTL <- setdiff(NR_respdeg, all_QTLSknow) #4293
unique(non_eQTL)

toplistall %>%
  # dplyr::filter(as.numeric(adj.P.Val)<0.1) %>%
  filter(ENTREZID %in% non_eQTL) %>%
  mutate(id = as.factor(id)) %>%
  mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%
  ggplot(., aes(x=id, y =abs(logFC)))+
  geom_boxplot(aes(fill=id))+
  fill_palette(palette =drug_palNoVeh)+
  guides(fill=guide_legend(title = "Treatment"))+
  #facet_wrap(~time,labeller = (time = facettimelabel) )+
  theme_bw()+
  xlab("")+
  ylab("Log Fold Change")+
  theme_bw()+
  facet_wrap(~time)+
  ggtitle("non-eQTL sigDEGs")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        # axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        strip.background = element_rect(fill = "transparent"),
        axis.text = element_text(size = 8, color = "black", angle = 0),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+


  geom_signif(comparisons = list(c("Daunorubicin", "Doxorubicin"),

                                 c("Epirubicin", "Doxorubicin"),
                                 c("Mitoxantrone", "Doxorubicin"),
                                c("Doxorubicin", "Trastuzumab")),

              test = "t.test",
              map_signif_level = FALSE,step_increase = 0.1,
              textsize = 4)

#length(intersect(backGL$ENTREZID,non_eQTL))
# toplistall %>%
#   filter(ENTREZID %in%knowles4) %>%
#   mutate(id = as.factor(id)) %>%
#   mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%
#   ggplot(., aes(x=id, y =abs(logFC)))+
#   geom_boxplot(aes(fill=id))+
#   fill_palette(palette =drug_palNoVeh)+
#   guides(fill=guide_legend(title = "Treatment"))+
#   #facet_wrap(~time,labeller = (time = facettimelabel) )+
#   theme_bw()+
#   xlab("")+
#   ylab(" |Log Fold Change|")+
#   theme_bw()+
#   facet_wrap(~time)+
#   ggtitle("non-eQTL n= 14044")+
#   theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
#         axis.title = element_text(size = 15, color = "black"),
#         # axis.ticks = element_line(linewidth = 1.5),
#         axis.line = element_line(linewidth = 1.5),
#         strip.background = element_rect(fill = "transparent"),
#         axis.text = element_text(size = 8, color = "black", angle = 0),
#         strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
#
#
#   geom_signif(comparisons = list(c("Daunorubicin", "Doxorubicin"),
#
#                                  c("Epirubicin", "Doxorubicin"),
#                                  c("Mitoxantrone", "Doxorubicin"),
#                                  c("Trastuzumab","Doxorubicin")),
#
#               test = "t.test",
#               map_signif_level = FALSE,step_increase = 0.1,
#               textsize = 4)
#


# doing just the overlapped genes -----------------------------------------
 sig_DEGname.noqtl <- setdiff(all_sigDEGname,all_QTLSknow)#8260
 all_QTLSknow45<- union(knowles5,knowles4)

 non_eQTL45 <- setdiff(all_sigDEGname, all_QTLSknow45)
 #8314
 toplistall %>%

   filter(ENTREZID %in% non_eQTL45) %>%
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
   ggtitle("Excluding K4 and K5 eQTLs")+
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

#other gene sets:  GWAS from counts cpm ----------------------------------------------------------
 #cpmcounts <- read.csv("data/filtered_cpm_counts.csv", row.names = 1)

 time <- rep((rep(c("3h", "24h"), c(6,6))), 6)
 time <- ordered(time, levels =c("3h", "24h"))
 drug <- rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Vehicle"),12)
 indv <- as.factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))
 labeltop <- (interaction(substring(drug, 0, 2), indv, time))

 #genesof interest:
 GWASsnps <- c('RARG', 'ZNF740', 'ITGB7',  'PRDM2', 'KAZN', 'GDF5', 'FRS2', 'HDDC2', 'EEF1B2', 'SLC28A3')
 totalsnpgenes <- getBM(attributes=my_attributes,filters ='hgnc_symbol',
                        values = GWASsnps, mart = ensembl)
replacename <- unique(totalsnpgenes[,c(1,3)])


 norm_GWAS <- cpmcounts[rownames(cpmcounts)%in% totalsnpgenes$entrezgene_id,]
 colnames(norm_GWAS) <- labeltop
 norm_GWAS %>%rownames_to_column("Gene") %>%
   pivot_longer(-Gene, names_to = "treatment",values_to = "counts") %>%
 separate_wider_delim(treatment,".",names =c("Drug", "indv", "time")) %>%
 # mutate(Drug =case_match(Drug, "Da"~"Daunorubicin",
 #                         "Do"~"Doxorubicin",
 #                         "Ep"~"Epirubicin",
 #                         "Mi"~"Mitoxantrone",
 #                         "Tr"~"Trastuzumab","Ve" ~"Vehicle", .default = Drug)) %>%
   mutate(time=factor(time, levels =c("3h", "24h"))) %>%
   mutate(indv=factor(indv)) %>%
   mutate(Gene=case_match(Gene, "1933"~"EEF1B2",
                          "10818" ~"FRS2",
                          "8200" ~ "GDF5",
          "51020"  ~     "HDDC2",
          "3695"  ~     "ITGB7",
          "23254"  ~      "KAZN",
           "7799"  ~     "PRDM2",
          "5916"   ~     "RARG",
         "64078"  ~   "SLC28A3",
         "283337"~"ZNF740")) %>%
 ggplot(., aes(x=Drug, y=counts))+
   geom_boxplot(position="identity",aes(fill=Drug))+
   geom_point(aes(col=indv, size=1.5, alpha=0.5))+
   guides(alpha= "none", size= "none")+
   scale_color_brewer(palette = "Dark2")+
   scale_fill_manual(values=drug_palc)+
   facet_wrap(Gene~time,  nrow=6, ncol=4)+
   theme_bw()+
   ylab(expression("Log "[2]~"cpm "))+
   xlab("")+
   xlab("")+
   theme(strip.background = element_rect(fill = "transparent"),
         plot.title = element_text(size=18,hjust = 0.5),
         axis.title = element_text(size = 15, color = "black"),
         axis.ticks = element_line(linewidth = 1.5),
         axis.line = element_line(linewidth = 1.5),
         axis.text.x = element_text(size = 12, color = "black", angle = 0),
         strip.text.x = element_text(size = 10, color = "black", face = "bold.italic"))
GwassigVDA24 <- sigVDA24%>%
  filter(ENTREZID %in% totalsnpgenes$entrezgene_id)
GwassigVDX24 <- sigVDX24%>%
  filter(ENTREZID %in% totalsnpgenes$entrezgene_id)
GwassigVEP24 <- sigVEP24%>%
  filter(ENTREZID %in% totalsnpgenes$entrezgene_id)
GwassigVMT24 <- sigVMT24%>%
  filter(ENTREZID %in% totalsnpgenes$entrezgene_id)

# Seoane 2019 gene list ---------------------------------------------------


 chrom_reg_Seoane <- read_csv(file = "data/Seonane2019supp1.txt",col_types = cols(...1 = col_skip()))

 Seoane_2019 <- chrom_reg_Seoane[,2]
 names(Seoane_2019) <- "ENTREZID"


 length(intersect(Seoane_2019$ENTREZID, NoResp$ENTREZID))


 norm_Seoane <- cpmcounts[rownames(cpmcounts)%in% Seoane_2019$ENTREZID,]
 colnames(norm_Seoane) <- labeltop

 Seoaneoverlap24 <- toplist24hours %>%
   dplyr::filter(ENTREZID%in%Seoane_2019$ENTREZID) %>%
   group_by(id) %>%
   filter(adj.P.Val<0.1) %>%
   count()

 NRseoane<- tibble_row( id = "No Response", n = 109)
 Seoaneoverlap <- merge(Seoaneoverlap,NRseoane, all =TRUE)[,1:2]
 Seoaneoverlap3 <- toplist3hours %>%
   dplyr::filter(ENTREZID%in%Seoane_2019$ENTREZID) %>%
   group_by(id) %>%
   filter(adj.P.Val<0.1) %>%
   count()



 ggplot(Seoaneoverlap3, aes(x=id, y=n))+
   geom_col(position = "dodge", aes(fill=id))+
   scale_color_brewer(palette = "Dark2",guide = "none")+
   scale_fill_manual(values=drug_palc[c(1:4,6)])+
   theme_bw()+

   ylab("Count")+
   xlab("")+
   ggtitle("Overlap with Seoane data")+
   theme(plot.title = element_text(size=18,hjust = 0.5),
         axis.title = element_text(size = 15, color = "black"),
         axis.ticks = element_line(linewidth = 1.5),
         axis.line = element_line(linewidth = 1.5),
         axis.text.x = element_text(size = 12, color = "black", angle = 0),
         strip.text.x = element_text(size = 15, color = "black", face = "bold"))











 norm_Seoane %>%rownames_to_column("Gene") %>%
   pivot_longer(-Gene, names_to = "treatment",values_to = "counts") %>%
   separate_wider_delim(treatment,".",names =c("Drug", "indv", "time")) %>%
   # mutate(Drug =case_match(Drug, "Da"~"Daunorubicin",
   #                         "Do"~"Doxorubicin",
   #                         "Ep"~"Epirubicin",
   #                         "Mi"~"Mitoxantrone",
   #                         "Tr"~"Trastuzumab","Ve" ~"Vehicle", .default = Drug)) %>%
   mutate(time=factor(time, levels =c("3h", "24h"))) %>%
   mutate(indv=factor(indv)) %>%
   filter(Drug=="Da"|Drug == "Ve") %>%

   # mutate(Gene=case_match(Gene, "1933"~"EEF1B2","10818" ~"FRS2",
   #                        "8200" ~       "GDF5",
   #                        "51020"  ~     "HDDC2",
   #                        "3096"  ~    "HIVEP1",
   #                        "3695"  ~     "ITGB7",
   #                        "23254"  ~      "KAZN",
   #                        "7799"  ~     "PRDM2",
   #                        "5916"   ~     "RARG",
   #                        "64078"  ~   "SLC28A3", .default = Gene)) %>%
   ggplot(., aes(x=Drug, y=counts))+
   geom_boxplot(position="identity",aes(fill=Drug))+
   geom_point(aes(col=indv, size=1.5, alpha=0.5))+
   guides(alpha= "none", size= "none")+
   scale_color_brewer(palette = "Dark2",guide = "none")+
   scale_fill_manual(values=drug_pal)+
   #facet_wrap_paginate(Gene~time, scales="free_y", nrow=6, ncol=6)+
   facet_wrap(~"time") +
   theme_bw()+
   ylab(expression("Log "[2]~"cpm "))+
   xlab("")+
   xlab("")+
   theme(strip.background = element_rect(fill = "transparent"),
         plot.title = element_text(size=18,hjust = 0.5),
         axis.title = element_text(size = 15, color = "black"),
         axis.ticks = element_line(linewidth = 1.5),
         axis.line = element_line(linewidth = 1.5),
         axis.text.x = element_text(size = 12, color = "black", angle = 0),
         strip.text.x = element_text(size = 10, color = "black", face = "bold"))


library(ggforce)

 for(i in 1:n_pages(p)){
   p_save <-  p +
     facet_wrap_paginate(~ Genes, ncol = 2, nrow = 2, page = i)
   ggsave(plot = p_save, filename = paste0('Downloads/page_', i, '.jpg'))
 }

# number of eQTLS by sigDA,etc at 24 hours --------------------------------

##24hours no response
 #first combine all 3 and 24 hour deg sets
 sig3all<- union(sigVMT3$ENTREZID,union(sigVEP3$ENTREZID,union(sigVDA3$ENTREZID,sigVDX3$ENTREZID)))
 sig24all<- union(sigVMT24$ENTREZID,union(sigVEP24$ENTREZID,union(sigVDA24$ENTREZID,sigVDX24$ENTREZID)))
 length(unique(sig24all))#9294
 sig_3_24_all <- union(sig3all,sig24all)  #9491
 ##now remove all expressed sig degs from backGL no make the NRresp list
 NR_respdeg <- setdiff(backGL$ENTREZID, sig_3_24_all)  #4593
length(intersect(NR_respdeg, knowles4)) #221
length(intersect(sigVDA24$ENTREZID,knowles4))#239
length(intersect(sigVDX24$ENTREZID,knowles4))#220
length(intersect(sigVEP24$ENTREZID,knowles4))#206
length(intersect(sigVMT24$ENTREZID,knowles4))#47

length(intersect(sigVDA24$ENTREZID,knowles4))#239
length(intersect(sigVDX24$ENTREZID,knowles4))#220
length(intersect(sigVEP24$ENTREZID,knowles4))#206
length(intersect(sigVMT24$ENTREZID,knowles4))#47



 non_eQTL45 <- setdiff(all_sigDEGname, all_QTLSknow45)

Overlapk4 <- toplistall %>%
  dplyr::filter(ENTREZID%in%knowles4) %>%
  group_by(id,time) %>%
  filter(adj.P.Val<0.1) %>%
 count()
Overlapk4 <- merge(Overlapk4,NR_respdeg, all =TRUE)[,1:2]
NR4<- tibble_row( id = "No Response", n = 248)

know4time <- tibble(id=c("No Response", "Daun","Doxo","Epi", "Mito", "Tras"),
                    time = c("3hour","24hour"),
                    n=c(6,0,6,2,0,331,217,200,)))

intersect(NoResp$ENTREZID,knowles4)

Overlapk5 <- toplist24hours %>%
  dplyr::filter(ENTREZID%in%knowles5) %>%
  group_by(id) %>%
  # filter(adj.P.Val<0.1) %>%
  count()
intersect(NoResp$ENTREZID,knowles5)
NRk5 <- tibble_row( id = "No Response", n = 141)
Overlapk5 <- merge(Overlapk5,NRk5, all =TRUE)[,1:2]


length(intersect(intersect(knowles4, knowles5), NoResp$ENTREZID))

drug_palc <- c("#8B006D","#DF707E","#F1B72B", "#3386DD","#707031","#41B333")

ggplot(Overlapk4, aes(x=id, y=n))+
  geom_col(aes(fill=id))+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_palc[c(1:4,6)])+
  theme_bw()+
  ylab("number of Overlaps with Knowles 4")+
  xlab("")+
  xlab("")+
  theme(plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 12, color = "white", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))



ggplot(Overlapk5, aes(x=id, y=n))+
  geom_col(aes(fill=id))+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_palc[c(1:4,6)])+
  theme_bw()+
  ylab("Count")+
  xlab("")+
  ggtitle("Overlaps with Knowles 5")+
  theme(plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 12, color = "black", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))


##make a set###---------------------------------------------------------

# this is the fill model of just the distribution of counts


knover45 <- list(Overlapk5, Overlapk4)
names(knover45) <- c("Knowles 4", "Knowles 5")
knover45 <- bind_rows(list("all eQTLs"=Overlapk4,"all response eQTLs"=Overlapk5), .id ="Overlap")
knover45$counts <- geneset_length <- rep(c(length(sigVDA24$ENTREZID),length(sigVDX24$ENTREZID),length(sigVEP24$ENTREZID), length(sigVMT24$ENTREZID),length(NoResp$ENTREZID)),2)


knover45 %>%
  # ggplot(., aes(counts, y=n))+
  ggplot(., aes(x=Overlap, y=n))+
  geom_col(position = "fill", aes(fill=id))+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_palc[c(1:4,6)])+
  theme_bw()+

  ylab("Count")+
  xlab("")+
  ggtitle("distribution of all QTLs and response eQTLS")+
  theme(plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 12, color = "black", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))



###proportion ------------------------------------------------------------
knover45 %>%
  mutate(percent =n/counts) %>%
  # ggplot(., aes(counts, y=n))+
  ggplot(., aes(x=Overlap, y=percent))+
  geom_col(position = "fill", aes(fill=id))+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_palc[c(1:4,6)])+
  theme_bw()+
    ylab("Percent of expressed genes")+
  xlab("")+
  ggtitle("Percent of reQTL distribution")+
  theme(plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 12, color = "black", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))



# plot the nonDE stuff ----------------------------------------------------






intersect(knowles5,toplist24hours$ENTREZID)
# https://bitbucket.org/jdblischak/tb/src/master/reQTL.Rmd

k45NRvenn <- list(as.numeric(knowles4), as.numeric(knowles5), NoResp$ENTREZID)

ggVennDiagram::ggVennDiagram(k45NRvenn,
              category.names = c("K4 n= 521","K5 n=377", "No response n=5776"),
              show_intersect = FALSE,
              set_color = "black",
              label = "both",
              label_percent_digit = 1,
              label_size = 3.5,
              label_alpha = 0,
              label_color = "black",
              edge_lty = "solid") +
  scale_x_continuous(expand = expansion(mult = .2))+
  scale_fill_gradientn(colours = terrain.colors(7))+
  #scale_fill_gradient(low = "red2", high = "yellow")+
  labs(title = "Knowles comparison to NR", caption = "")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))
list3totvenn <- get.venn.partitions(total3)



# boxplot absolute Fold change --------------------------------------------

#load lists needed:  Knowles4, Knowles5, and all expressed genes without the 4 and 5 sets (take out the all_QTLSknow)
cpmcounts <- read.csv("data/filtered_cpm_counts.csv", row.names = 1)


toplistall %>%
  mutate(id = as.factor(id)) %>%
  mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%
  filter(ENTREZID %in% knowles4) %>%
  ggplot(., aes(x=id, y=abs(logFC)))+
  geom_boxplot(aes(fill=id))+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_palc[c(1:4,6)])+
  theme_bw()+

  ylab("|log FC|")+
  xlab("")+
  facet_wrap(~time)+
  ggtitle("abs log FC of all eQTLs in data set")+
  theme(plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 12, color = "black", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))

toplistall %>%
  filter(ENTREZID %in% knowles5) %>%
  ggplot(., aes(x=id, y=abs(logFC)))+
  geom_boxplot(aes(fill=id))+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_palc[c(1:4,6)])+
  theme_bw()+

  ylab("|log FC|")+
  xlab("")+
  facet_wrap(~time)+
  ggtitle("abs log FC of all reQTL across treatments")+
  theme(plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 12, color = "black", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))

toplistall %>%
  mutate(id = as.factor(id)) %>%
  mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%

  filter(ENTREZID %in% motif_NR) %>%
  ggplot(., aes(x=id, y=abs(logFC)))+
  geom_boxplot(aes(fill=id))+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_palc[c(1:4,6)])+
  theme_bw()+
  ylab("|log FC|")+
  xlab("")+
  facet_wrap(~time)+
  ggtitle("abs log FC no-response genes motif1")+
  theme(plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 12, color = "black", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))
  geom_signif(comparisons = list(c("Daunorubicin", "Doxorubicin"),
                                 c("Epirubicin", "Doxorubicin"),
                                 c("Mitoxantrone", "Doxorubicin"),
                                 c("Trastuzumab","Doxorubicin")),
              test = "t.test",
              map_signif_level = FALSE,step_increase = 0.1,
              textsize = 4)

  toplistall %>%
    mutate(id = as.factor(id)) %>%
    mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%

    filter(ENTREZID %in% motif3_TI) %>%
    ggplot(., aes(x=id, y=abs(logFC)))+
    geom_boxplot(aes(fill=id))+
    scale_color_brewer(palette = "Dark2",guide = "none")+
    scale_fill_manual(values=drug_palc[c(1:4,6)])+
    theme_bw()+
    ylab("|log FC|")+
    xlab("")+
    facet_wrap(~time)+
    ggtitle("abs log FC time-independent")+
    theme(plot.title = element_text(size=18,hjust = 0.5),
          axis.title = element_text(size = 15, color = "black"),
          axis.ticks = element_line(linewidth = 1.5),
          axis.line = element_line(linewidth = 1.5),
          axis.text.x = element_text(size = 12, color = "black", angle = 0),
          strip.text.x = element_text(size = 15, color = "black", face = "bold"))


  toplistall %>%
    mutate(id = as.factor(id)) %>%
    mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%

    filter(ENTREZID %in% motif4_LR) %>%
    ggplot(., aes(x=id, y=abs(logFC)))+
    geom_boxplot(aes(fill=id))+
    scale_color_brewer(palette = "Dark2",guide = "none")+
    scale_fill_manual(values=drug_palc[c(1:4,6)])+
    theme_bw()+
    ylab("|log FC|")+
    xlab("")+
    facet_wrap(~time)+
    ggtitle("abs log FC Late response")+
    theme(plot.title = element_text(size=18,hjust = 0.5),
          axis.title = element_text(size = 15, color = "black"),
          axis.ticks = element_line(linewidth = 1.5),
          axis.line = element_line(linewidth = 1.5),
          axis.text.x = element_text(size = 12, color = "black", angle = 0),
          strip.text.x = element_text(size = 15, color = "black", face = "bold"))

  toplistall %>%
    mutate(id = as.factor(id)) %>%
    mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%

    filter(ENTREZID %in% motif5_ER) %>%
    ggplot(., aes(x=id, y=abs(logFC)))+
    geom_boxplot(aes(fill=id))+
    scale_color_brewer(palette = "Dark2",guide = "none")+
    scale_fill_manual(values=drug_palc[c(1:4,6)])+
    theme_bw()+
    ylab("|log FC|")+
    xlab("")+
    facet_wrap(~time)+
    ggtitle("abs log FC Early response")+
    theme(plot.title = element_text(size=18,hjust = 0.5),
          axis.title = element_text(size = 15, color = "black"),
          axis.ticks = element_line(linewidth = 1.5),
          axis.line = element_line(linewidth = 1.5),
          axis.text.x = element_text(size = 12, color = "black", angle = 0),
          strip.text.x = element_text(size = 15, color = "black", face = "bold"))




# cormotif distrubution ---------------------------------------------------
DEG_cormotif <- readRDS("data/DEG_cormotif.RDS")
motif_NR <- DEG_cormotif$motif_NR
motif_TI <- DEG_cormotif$motif_TI
motif_LR <- DEG_cormotif$motif_LR
motif_ER <- DEG_cormotif$motif_ER

length(intersect(knowles5, motif_NR))
length(intersect(knowles5, motif_TI))
length(intersect(knowles5, motif_LR))
length(intersect(knowles5, motif_ER))

Overlapk4 <- merge(Overlapk1,NR4, all =TRUE)[,1:2]

cormotfi4<- tibble( id =c("No response","Time-Independent response", "Late response","Early response"), n = c(329, 3,138,4), totalgene= c(7362,432,4850,length(motif_ER)))
cormotfi5 <- tibble( id =c("No response","Time-Independent response", "Late response","Early response"), n = c(199,12,139,9), totalgene= c(length(motif_NR),length(motif_TI),length(motif_LR), length(motif_ER)))



cormotfi5 %>%
  mutate(percent = n/totalgene) %>%
  mutate(id=factor(id, levels = c("Time-Independent response",
                                  "Early response","Late response","No response"))) %>%
ggplot(., aes(x=id, y=n))+
  geom_col(aes(fill=id))+
  scale_fill_brewer(guide = "none")+

  theme_bw()+
  ylab("Count")+
  xlab("")+
  ggtitle("reQTL distribution in cormotif set")+
  theme(plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 12, color = "black", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))



cormotfi4 %>%
  mutate(percent = n/totalgene) %>%
  mutate(id=factor(id, levels = c("Time-Independent response",
                                  "Early response","Late response","No response"))) %>%
  ggplot(., aes(x=id, y=percent))+
  geom_col()+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_palc[c(1:4,6)])+
  theme_bw()+
  ylab("Percent of gene set")+
  xlab("")+
  ggtitle("marginal eQTLs between sets")+
  theme(plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 12, color = "black", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))








#
# Overlap3 <- toplist24hours %>%
#   dplyr::filter(ENTREZID %in% knowles4) %>%
#   dplyr::filter(ENTREZID %in% motif3_TI) %>%
#   group_by(id) %>%
#     count()
# Overlap4 <- toplist24hours %>%
#   dplyr::filter(ENTREZID %in% knowles4) %>%
#   dplyr::filter(ENTREZID %in% motif4_LR) %>%
#   group_by(id) %>%
#    count()
# Overlap5 <- toplist24hours %>%
#   dplyr::filter(ENTREZID%in% motif5_ER) %>%
#   group_by(id) %>%
#   filter(adj.P.Val<0.1) %>%
#   count()
# cor_qtl <-bind_rows(list(Overlap1, Overlap3, Overlap4, Overlap5), .id = "motif")#, all = TRUE)[,1:2]
#
# backGL <- read.csv("data/backGL.txt")
# NRresp <- read_csv("data/cormotif_NRset.txt")
##first filter and
