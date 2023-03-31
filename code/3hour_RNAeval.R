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
drug_palNoVeh <- c("#8B006D","#DF707E","#F1B72B", "#3386DD","#707031","#41B333")
#named colors: dark pink,Red,yellow,blue, dark grey, green
# yarrr::piratepal("all")
# indv_pal <- c(appletv color pallett from yarrr package or southpark)
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
  ylim(-5,5)+
  xlim(-5,5)+
  stat_cor(aes(label = after_stat(rr.label)),
           color = "red",
           geom = "label",
           label.y = 4, label.x=-4)+
  ggtitle("Trastuzumab and Doxorubicin 3 hours")

toplist3hours %>%
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
  ggtitle("Mitoxantrone and Doxorubicin 3 hours")


toplist3hours %>%
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
  ggtitle("Epirubicin and Doxorubicin 3 hours")

toplist3hours %>%
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
           geom = "label",
           label.y = 4, label.x=-4)+
  ggtitle("Daunorubicin and Doxorubicin 3 hours")

df <- data.frame(x = 1:10, y = 1)


my_barcolor<- c("#FFFFFF" ,"#E6E6E6", "#CCCCCC" ,"#B3B3B3" ,"#999999" ,"#808080", "#666666",
"#4C4C4C" ,"#333333" ,"#191919" ,"#000000")
 ggplot(df, aes(x, y)) +
  geom_tile(aes(fill = x)) +
  scale_fill_gradientn(colors = my_barcolor, guide = "colorbar")
