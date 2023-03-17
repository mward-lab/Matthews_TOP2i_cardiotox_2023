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


#Detach all  packages
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base","package:workflowr")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)

}

detachAllPackages()
####plotting the DXoggplots
efit2<- readRDS(file ="data/efit2results.RDS")

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


toplist24hours %>%
  filter(id =="Trastuzumab"|id=="Doxorubicin") %>%
  dplyr::select(id,logFC,ENTREZID) %>%
  pivot_wider(names_from =id, values_from = logFC) %>%
  ggplot(.,aes(x= Doxorubicin, y=Trastuzumab))+
  geom_point()+
  geom_smooth(method="lm")+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()+
  stat_cor(aes(label = after_stat(rr.label)),
           color = "red",
           geom = "label",
           label.y = 5, label.x=-8)+
  ggtitle("TR and DX 24")

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



  plot_pathway(
    data = df,
    comp.names = NULL,
    gene.id.type = "ENSEMBL",
    FC.cutoff = 1.3,
    FDR.cutoff = 0.05,
    FCflag = "logFC",
    FDRflag = "adj.P.Val",
    Fisher.cutoff = 0.1,
    Fisher.up.cutoff = 0.1,
    Fisher.down.cutoff = 0.1,
    plot.save.to = NULL,
    pathway.db = "rWikiPathways"
  )
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









