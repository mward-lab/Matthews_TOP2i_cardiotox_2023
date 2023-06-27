#go analysis of stuffs
library(readxl)
library(rstatix)
library(tidyverse)
library(zoo)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(dplyr)
library(car)
library(limma)

drug_palNoVeh <- c("#8B006D","#DF707E","#F1B72B", "#3386DD","#707031","#41B333")
#named colors: dark pink,Red,yellow,blue, dark grey, green
#Detach all  packages
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base","package:workflowr")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)

}

detachAllPackages()
####plotting the DXoggplots
efit2<- readRDS(file ="data/efit2_final.RDS")

# saveRDS(efit2, file ="data/efit2results.RDS")


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
saveRDS(toplistall, "data/toplistall.RDS")
readRDS("data/toplistall.RDS")
toplist24hours %>%
  filter(id =="Daunorubicin"|id=="Doxorubicin") %>%
  dplyr::select(id,logFC,ENTREZID) %>%
  pivot_wider(names_from =id, values_from = logFC) %>%
  ggplot(.,aes(x= Doxorubicin, y=Daunorubicin))+
  geom_point()+
  geom_smooth(method="lm")+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()+
  ylim(-5,5)+
  xlim(-5,5)+
  stat_cor(aes(label = after_stat(rr.label)),
           color = "red",
           geom = "text",
           label.y = 4, label.x=-4)+
  ggtitle("Daunorubicin and Doxorubicin 24 hours")



toplist24hours %>%
  filter(id =="Epirubicin"|id=="Doxorubicin") %>%
  dplyr::select(id,logFC,ENTREZID) %>%
  pivot_wider(names_from =id, values_from = logFC) %>%
  ggplot(.,aes(x= Doxorubicin, y=Epirubicin))+
  geom_point()+
  geom_smooth(method="lm")+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()+
  ylim(-5,5)+
  xlim(-5,5)+
  stat_cor(aes(label = after_stat(rr.label)),
           color = "red",
           geom = "label",
           label.y = 4, label.x=-4)+
  ggtitle("Epirubicin and Doxorubicin 24 hours")


toplist24hours %>%
  filter(id =="Mitoxantrone"|id=="Doxorubicin") %>%
  dplyr::select(id,logFC,ENTREZID) %>%
  pivot_wider(names_from =id, values_from = logFC) %>%
  ggplot(.,aes(x= Doxorubicin, y=Mitoxantrone))+
  geom_point()+
  geom_smooth(method="lm")+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()+
  ylim(-5,5)+
  xlim(-5,5)+
  stat_cor(aes(label = after_stat(rr.label)),
           color = "red",
           geom = "label",
           label.y = 4, label.x=-4)+
  ggtitle("Mitoxantrone and Doxorubicin 24 hours")




toplist24hours %>%
  filter(id =="Trastuzumab"|id=="Doxorubicin") %>%
  dplyr::select(id,logFC,ENTREZID) %>%
  pivot_wider(names_from =id, values_from = logFC) %>%
  ggplot(.,aes(x= Doxorubicin, y=Trastuzumab))+
  geom_point()+
  geom_smooth(method="lm")+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()+
  ylim(-5,5)+
  xlim(-5,5)+
  stat_cor(aes(label = after_stat(rr.label)),
           color = "red",
           geom = "label",
           label.y = 4, label.x=-4)+
  ggtitle("Trastuzumab and Doxorubicin 24 hours")





sigVDX24 %>% inner_join(sigVDA24, by="ENTREZID") %>%
  ggplot(., aes(x=logFC.x, y=logFC.y, group=ENTREZID))+
  geom_point()+
  geom_smooth(method="lm")+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()+
  stat_cor(aes(label = after_stat(rr.label)),
           color = "red",
           geom = "label",
           label.y = 5, label.x=-8)+
  ggtitle("DA and DX 24")

sigVDX24 %>% inner_join(sigVEP24, by="ENTREZID") %>%
  ggplot(., aes(x=logFC.x, y=logFC.y, group=ENTREZID))+
  geom_point()+
  geom_smooth(method="lm")+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()+
  stat_cor(aes(label = after_stat(rr.label)),
           color = "red",
           geom = "label",
           label.y = 5, label.x=-8)+
  ggtitle("EP and DX 24")
sigVDX24 %>% inner_join(sigVMT24, by="ENTREZID") %>%
  ggplot(., aes(x=logFC.x, y=logFC.y, group=ENTREZID))+
  geom_point(aes(alpha=0.5))+
  geom_smooth(method="lm")+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()+
  stat_cor(aes(label = after_stat(rr.label)),
           color = "red",
           geom = "label",
           label.y = 5, label.x=-8)+
  ggtitle("MT and DX 24")

#
#
#   plot_pathway(
#     data = df,
#     comp.names = NULL,
#     gene.id.type = "ENSEMBL",
#     FC.cutoff = 1.3,
#     FDR.cutoff = 0.05,
#     FCflag = "logFC",
#     FDRflag = "adj.P.Val",
#     Fisher.cutoff = 0.1,
#     Fisher.up.cutoff = 0.1,
#     Fisher.down.cutoff = 0.1,
#     plot.save.to = NULL,
#     pathway.db = "rWikiPathways"
#   )
  toplistall %>%
    filter(id =="Trastuzumab"|id=="Doxorubicin") %>%
    dplyr::select(time,id,logFC,ENTREZID) %>%
    pivot_wider(names_from =id, values_from = logFC) %>%
    ggplot(.,aes(x= Doxorubicin, y=Trastuzumab,group=time))+
    geom_point()+
    geom_smooth(method="lm")+
    geom_abline(slope = 1, intercept = 0)+
   facet_wrap('time',ncol=1,nrow=)+
    theme_bw()+
    stat_cor(aes(label = after_stat(rr.label)),
             color = "red",
             geom = "label",
             label.y = 5, label.x=-8)+
    ggtitle("TR and DX by time")


# working on FC  ----------------------------------------------------------

toplist24hours %>% dplyr::filter(P.Value<0.1) %>%
    ggplot(.,aes(x=id, y = logFC))+
    geom_boxplot(position = "identity")+
    geom_signif(comparisons = list(c("Daunorubicin", "Doxorubicin"),
                                   c("Epirubicin", "Doxorubicin"),
                                   c("Mitoxantrone", "Doxorubicin"),
                                   c("Trastuzumab","Doxorubicin")),
                test = "t.test",
                map_signif_level = TRUE,
                step_increase = 0.1,
                textsize = 4)

  toplistall %>%
    filter(ENTREZID %in% DDEMresp) %>%
    mutate(id = as.factor(id)) %>%
    mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%
    ggplot(., aes(x=id, y =logFC))+
    geom_boxplot(aes(fill=id))+
    fill_palette(palette =drug_palNoVeh)+
    guides(fill=guide_legend(title = "Treatment"))+
    facet_wrap(~time,labeller = (time = facettimelabel) )+
    theme_bw()+
    xlab("")+
    ylab("Log Fold Change")+
    theme_bw()+
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
                map_signif_level = TRUE,step_increase = 0.1,
                textsize = 4)





##modified for poster


time.label <- c("3 hours","24 hours")
names(time.label) <- c('3_hours','24_hours')
  motif3_TI
  toplistall %>%

    mutate(id = as.factor(id)) %>%
    mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%
    ggplot(., aes(x=id, y =logFC))+
    # geom_hline(yintercept = 0,lty =2, linewidth=2, color= "red") +
    geom_boxplot(aes(fill=id))+
    fill_palette(palette =drug_palNoVeh)+
    guides(fill=guide_legend(title = "Treatment"))+

    facet_wrap(~time)+
    theme_bw()+
    xlab("")+
    ylab("Log Fold Change")+
    theme_bw()+
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
                map_signif_level = TRUE,step_increase = 0.05,
                textsize = 4)


  toplistall %>%
    filter(ENTREZID %in% motif5_ER) %>%
    mutate(id = as.factor(id)) %>%
    mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%
    ggplot(., aes(x=id, y =logFC))+
    geom_boxplot(aes(fill=id))+
    fill_palette(palette =drug_palNoVeh)+
    guides(fill=guide_legend(title = "Treatment"))+

    facet_wrap(~time, )+
    theme_bw()+
    xlab("Early Response motif 5")+
    ylab("Log Fold Change")+
    theme_bw()+
    theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = 15, color = "black"),
          # axis.ticks = element_line(linewidth = 1.5),
          axis.line = element_line(linewidth = 1.5),
          strip.background = element_rect(fill = "transparent"),
          axis.text = element_text(size = 8, color = "white", angle = 0),
          strip.text.x = element_text(size = 12, color = "black", face = "bold"))+


    geom_signif(comparisons = list(c("Daunorubicin", "Doxorubicin"),

                                   c("Epirubicin", "Doxorubicin"),
                                   c("Mitoxantrone", "Doxorubicin"),
                                   c("Trastuzumab","Doxorubicin")),
                test = "t.test",
                map_signif_level = TRUE,step_increase = 0.1,
                textsize = 4)
  toplistall %>%
    filter(ENTREZID %in% motif4_LR) %>%
    mutate(id = as.factor(id)) %>%
    mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%
    ggplot(., aes(x=id, y =logFC))+
    geom_boxplot(aes(fill=id))+
    fill_palette(palette =drug_palNoVeh)+
    guides(fill=guide_legend(title = "Treatment"))+

    facet_wrap(~time, )+
    theme_bw()+
    xlab("Late Response motif 4")+
    ylab("Log Fold Change")+
    theme_bw()+
    theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = 15, color = "black"),
          # axis.ticks = element_line(linewidth = 1.5),
          axis.line = element_line(linewidth = 1.5),
          strip.background = element_rect(fill = "transparent"),
          axis.text = element_text(size = 8, color ="white", angle = 0),
          strip.text.x = element_text(size = 12, color = "black", face = "bold"))+


    geom_signif(comparisons = list(c("Daunorubicin", "Doxorubicin"),

                                   c("Epirubicin", "Doxorubicin"),
                                   c("Mitoxantrone", "Doxorubicin"),
                                   c("Trastuzumab","Doxorubicin")),
                test = "t.test",
                map_signif_level = TRUE,step_increase = 0.1,
                textsize = 4)

  toplistall %>%
    filter(ENTREZID %in% motif1_NR) %>%
    mutate(id = as.factor(id)) %>%
    mutate(time=factor(time, levels=c("3_hours","24_hours"))) %>%
    ggplot(., aes(x=id, y =logFC))+
    geom_boxplot(aes(fill=id))+
    fill_palette(palette =drug_palNoVeh)+
    guides(fill=guide_legend(title = "Treatment"))+

    facet_wrap(~time, )+
    theme_bw()+
    xlab("No Response motif 1")+
    ylab("Log Fold Change")+
    theme_bw()+
    theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = 15, color = "black"),
          # axis.ticks = element_line(linewidth = 1.5),
          axis.line = element_line(linewidth = 1.5),
          strip.background = element_rect(fill = "transparent"),
          axis.text = element_text(size = 8, color = "white", angle = 0),
          strip.text.x = element_text(size = 12, color = "black", face = "bold"))+


    geom_signif(comparisons = list(c("Daunorubicin", "Doxorubicin"),

                                   c("Epirubicin", "Doxorubicin"),
                                   c("Mitoxantrone", "Doxorubicin"),
                                   c("Trastuzumab","Doxorubicin")),
                test = "t.test",
                map_signif_level = TRUE,step_increase = 0.1,
                textsize = 4)

  ##using DDEMresp list n=952
# to get gene graph -------------------------------------------------------
#pick the list you want, reccommended to limit to 5 genes per list
  # wantedgl <-








  toplistall %>% filter(adj.P.Val < .1) %>%

    filter(id =="Epirubicin"|id=="Doxorubicin") %>%
      dplyr::select(time,id,logFC,ENTREZID) %>%
      pivot_wider(names_from =id, values_from = logFC) %>%
      ggplot(.,aes(x= Doxorubicin, y=Epirubicin,group=time))+
      geom_point()+
      geom_smooth(method="lm")+
      geom_abline(slope = 1, intercept = 0)+
      facet_wrap('time')+
      theme_bw()+
      stat_cor(aes(label = after_stat(rr.label)),
               color = "red",
               geom = "label")+
      ggtitle("Ep and DX by time")

  toplistall %>% filter(adj.P.Val < .1) %>%

    filter(id =="Daunorubicin"|id=="Doxorubicin") %>%
    dplyr::select(time,id,logFC,ENTREZID) %>%
    pivot_wider(names_from =id, values_from = logFC) %>%
    ggplot(.,aes(x= Doxorubicin, y=Daunorubicin,group=time))+
    geom_point()+
    geom_smooth(method="lm")+
    geom_abline(slope = 1, intercept = 0)+
    facet_wrap('time')+
    theme_bw()+
    stat_cor(aes(label = after_stat(rr.label)),
             color = "red",
             geom = "label")+
    ggtitle("DA and DX by time")

  toplistall %>%

    filter(id =="Daunorubicin"|id=="Doxorubicin") %>%
    dplyr::select(time,id,logFC,ENTREZID) %>%
    pivot_wider(names_from =id, values_from = logFC) %>%
    ggplot(.,aes(x= Doxorubicin, y=Daunorubicin,group=time))+
    geom_point()+
    geom_smooth(method="lm")+
    geom_abline(slope = 1, intercept = 0)+
    facet_wrap('time')+
    theme_bw()+
    stat_cor(aes(label = after_stat(rr.label)),
             color = "red",
             geom = "label")+
    ggtitle("DA and DX by time- no filter on adj p value")

  toplistall %>%
    filter(id =="Daunorubicin"|id=="Doxorubicin") %>%
    dplyr::select(time,id,logFC,ENTREZID) %>%
    pivot_wider(names_from =id, values_from = logFC) %>%
    ggplot(.,aes(x= Doxorubicin, y=Daunorubicin,group=time))+
    geom_point()+
    geom_smooth(method="lm")+
    geom_abline(slope = 1, intercept = 0)+
    facet_wrap('time')+
    theme_bw()+
    stat_cor(aes(label = after_stat(rr.label)),
             color = "red",
             geom = "label")+
    ggtitle("DA and DX by time- no filter on adj p value")

# GO of cormotifDEG -------------------------------------------------------


library(gprofiler2)

  DEG_cormotif <- readRDS("data/DEG_cormotif.RDS")
  motif1_NR <- DEG_cormotif$motif1_NR
  motif3_TI <- DEG_cormotif$motif3_TI
  motif4_LR <- DEG_cormotif$motif4_LR
  motif5_ER <- DEG_cormotif$motif5_ER




# TI gene set -------------------------------------------------------------



  gostresTI <- gost(query =motif3_TI,
                    organism = "hsapiens",
                    ordered_query = FALSE,
                    domain_scope = "custom",
                    measure_underrepresentation = FALSE,
                    evcodes = TRUE,
                    user_threshold = 0.01,
                    correction_method = c("fdr"),
                    custom_bg = backGL$ENTREZID,
                    sources="GO:BP", significant = FALSE)


  TI_gostresults <- gostplot(gostresTI, capped = FALSE, interactive = TRUE)

  tableTI <- gostresTI$result %>%
        dplyr::select(
      c(source, term_id, term_name,intersection_size, term_size, p_value)) #%>%
       # mutate_at(.vars = 6, .funs= scientific_format())  #use this  for table display only


  tableTI


  tableTI %>% dplyr::select(p_value,term_name,intersection_size) %>%
    slice_min(., n=20 ,order_by=p_value) %>%
    mutate(log_val = -log10(p_value)) %>%
   # slice_max(., n=10,order_by = p_value) %>%
   ggplot(., aes(x = log_val, y =reorder(term_name,p_value))) +
    geom_point(aes(size = intersection_size)) +
    ggtitle('Top2Bi-Time Independent enriched GO:BP terms') +
    xlab("-log 10 (p-value)")+
    ylab("GO: BP term")+
      theme_bw()



# ER_genes ----------------------------------------------------------------

  gostresER <- gost(query =motif5_ER,
                    organism = "hsapiens",
                    ordered_query = FALSE,
                    domain_scope = "custom",
                    measure_underrepresentation = FALSE,
                    evcodes = TRUE,
                    user_threshold = 0.01,
                    correction_method = c("fdr"),
                    custom_bg = backGL$ENTREZID,
                    sources="GO:BP", significant = FALSE)


  ER_gostresults <- gostplot(gostresER, capped = FALSE, interactive = TRUE)

  tableER <- gostresER$result %>%
    dplyr::select(
      c(source, term_id, term_name,intersection_size, term_size, p_value)) #%>%
  # mutate_at(.vars = 6, .funs= scientific_format())  #use this  for table display only


  tableER


  tableER %>% dplyr::select(p_value,term_name,intersection_size) %>%
    slice_min(., n=20 ,order_by=p_value) %>%
    mutate(log_val = -log10(p_value)) %>%
    # slice_max(., n=10,order_by = p_value) %>%
    ggplot(., aes(x = log_val, y =reorder(term_name,p_value))) +
    geom_point(aes(size = intersection_size)) +
    ggtitle('Top2Bi-Early response set enriched GO:BP terms') +
    xlab("-log 10 (p-value)")+
    ylab("GO: BP term")+
    theme_bw()


# LR genes ----------------------------------------------------------------

  gostresLR <- gost(query =motif4_LR,
                    organism = "hsapiens",
                    ordered_query = FALSE,
                    domain_scope = "custom",
                    measure_underrepresentation = FALSE,
                    evcodes = TRUE,
                    user_threshold = 0.01,
                    correction_method = c("fdr"),
                    custom_bg = backGL$ENTREZID,
                    sources="GO:BP", significant = FALSE)


  LR_gostresults <- gostplot(gostresLR, capped = FALSE, interactive = TRUE)

  tableLR <- gostresLR$result %>%
    dplyr::select(
      c(source, term_id, term_name,intersection_size, term_size, p_value)) #%>%
  # mutate_at(.vars = 6, .funs= scientific_format())  #use this  for table display only


  tableLR


  tableLR %>% dplyr::select(p_value,term_name,intersection_size) %>%
    slice_min(., n=20 ,order_by=p_value) %>%
    mutate(log_val = -log10(p_value)) %>%
    # slice_max(., n=10,order_by = p_value) %>%
    ggplot(., aes(x = log_val, y =reorder(term_name,p_value))) +
    geom_point(aes(size = intersection_size)) +
    ggtitle('Top2Bi-Late response set enriched GO:BP terms') +
    xlab("-log 10 (p-value)")+
    ylab("GO: BP term")+
    theme_bw()




# NR gene set -------------------------------------------------------------


  gostresNR <- gost(query =motif4_NR,
                    organism = "hsapiens",
                    ordered_query = FALSE,
                    domain_scope = "custom",
                    measure_underrepresentation = FALSE,
                    evcodes = TRUE,
                    user_threshold = 0.01,
                    correction_method = c("fdr"),
                    custom_bg = backGL$ENTREZID,
                    sources="GO:BP", significant = FALSE)


  NR_gostresults <- gostplot(gostresNR, capped = FALSE, interactive = TRUE)

  tableNR <- gostresNR$result %>%
    dplyr::select(
      c(source, term_id, term_name,intersection_size, term_size, p_value)) #%>%
  # mutate_at(.vars = 6, .funs= scientific_format())  #use this  for table display only


  tableNR


  tableNR %>% dplyr::select(p_value,term_name,intersection_size) %>%
    slice_min(., n=20 ,order_by=p_value) %>%
    mutate(log_val = -log10(p_value)) %>%
    # slice_max(., n=10,order_by = p_value) %>%
    ggplot(., aes(x = log_val, y =reorder(term_name,p_value))) +
    geom_point(aes(size = intersection_size, col=intersection_size)) +
    ggtitle('Top2Bi-No response set enriched GO:BP terms') +
    xlab("-log 10 (p-value)")+
    ylab("GO: BP term")+
    theme_bw()

  norm_counts <- read.csv("data/norm_counts.csv", row.names = 1)
norm_counts <- as.data.frame(log(cpm(norm_counts)))
colnames(norm_counts) <- labeltop

# graph by gene count -----------------------------------------------------

time <- rep((rep(c("3h", "24h"), c(6,6))), 6)
time <- ordered(time, levels =c("3h", "24h"))
drug <- rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Control"),12)
indv <- as.factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))
labelgenes sample

DEG_cormotif <- readRDS("data/DEG_cormotif.RDS")
motif1_NR <- DEG_cormotif$motif1_NR
motif3_TI <- DEG_cormotif$motif3_TI
motif4_LR <- DEG_cormotif$motif4_LR
motif5_ER <- DEG_cormotif$motif5_ER
# norm_counts <- read.csv("data/norm_counts.csv", row.names = 1)
# norm_counts <- as.data.frame(log(cpm(norm_counts)))
# colnames(norm_counts) <- labeltop
#


norm_counts %>% dplyr::filter(rownames(.)%in% motif3_TI) %>%
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





# GWASsnps ----------------------------------------------------------------

library(biomaRt)

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
my_chr <- c(1:22, 'M', 'X', 'Y')  ## creates a filter for each database
#attributes <- listAttributes((ensembl))
my_attributes <- c('entrezgene_id', 'ensembl_gene_id', 'hgnc_symbol')

GWASsnps <- c('RARG', 'HIVEP1', 'ITGB7',  'PRDM2', 'KAZN', 'GDF5', 'FRS2', 'HDDC2', 'EEF1B2', 'SLC28A3')
GWASsnps <- getBM(attributes=my_attributes,filters ='hgnc_symbol',
                    values = GWASsnps, mart = ensembl)
GWASsnps <- GWASsnps[unique(GWASsnps$entrezgene_id),]

toplist24hours %>%

