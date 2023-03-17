library(readxl)
library(rstatix)
library(tidyverse)
library(zoo)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(dplyr)
library(purrr)

library(edgeR)
###discrete colour scales!
scale_color_discrete()
drug_pal <- c("#8B006D","#DF707E","#F1B72B", "#3386DD","#707031","#41B333")
#named colors: dark pink,Red,yellow,blue, dark grey, green
yarrr::piratepal("all")
indv_pal <- c(appletv color pallett from yarrr package or southpark)
#Detach all  packages
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base","package:workflowr")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)

}

detachAllPackages()

# uploading and making 3 hour lists ---------------------------------------

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

toplist3hours %>%
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
  ggtitle("TR and DX 3")



# ggplot(RNAnormlist, aes(x=Drug, y=tnni))+
#   geom_boxplot(position = "identity")+
#   geom_point(aes(col=indv, size=2))+
#   geom_signif(comparisons = list(c("Daunorubicin", "Control"),
#                                  c("Doxorubicin", "Control"),
#                                  c("Epirubicin", "Control"),
#                                  c("Mitoxantrone", "Control"),
#                                  c("Trastuzumab","Control")),
#               test = "t.test",
#               map_signif_level = FALSE,step_increase = 0.1,
#               textsize = 6)+
#   ggtitle("Relative troponin I levels released in media")+
#   theme_bw()+
#   theme(strip.background = element_rect(fill = "transparent")) +
#   theme(
#     axis.title = element_text(size = 15, color = "black"),
#     axis.ticks = element_line(size = 1.5),
#     axis.text = element_text(size = 9, color = "black", angle = 0),
#     strip.text.x = element_text(size = 15, color = "black", face = "bold"))
#
#
