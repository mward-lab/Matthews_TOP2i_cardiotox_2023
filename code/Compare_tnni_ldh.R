
###script for comparing LDH and TNNI at 24 hours, and LDH at 24 and 48 hours at 0.5 uM
library(readxl)
library(rstatix)
library(tidyverse)
library(zoo)
library(ggplot2)
library(ggpubr)
library(ggsignif)
###discrete colour scales!
scale_color_brewer()
#indv   scale_color_brewer(palette = "Dark2")
drug_pal <- c("#8B006D","#DF707E","#F1B72B", "#3386DD","#707031","#41B333")
#

#using infofrom wilkelab.org  https://wilkelab.org/SDS375/slides/color-scales.html#47
  # scale_color_manual(
  #   values = c(West = "#E69F00", South = "#56B4E9", Midwest = "#009E73", Northeast = "#F0E442")
  # )

scale_color_brewer(palette = "Dark2")  # for use in indivduals
#named colors: dark pink,Red,yellow,blue, dark grey, green
yarrr::piratepal("all")
level_order <- c('71','75','77','78','79','87')
level_order <- c('71','75','77','78','79','87')
level_order2 <- c('75','87','77','79','78','71')
#southpark <- c(71 = "#2F86FFFF", 75 ="#EBAB16FF",77= "#DE0012FF", 78 ="#22C408FF", 79 ="#FECDAAFF", 87= "#F14809FF")
#Detach all  packages
# okabe ito colors
# orange	#E69F00	230, 159, 0
# sky blue	#56B4E9	86, 180, 233
# bluish green	#009E73	0, 158, 115
# yellow	#F0E442	240, 228, 66
# blue	#0072B2	0, 114, 178
# vermilion	#D55E00	213, 94, 0
# reddish purple	#CC79A7	204, 121, 167
# black	#000000	0, 0, 0

list(ggplot2.discrete.fill = okabe)




# TNNI 24 hour ------------------------------------------------------------


TNNIelisadata <- read_excel("~/Ward Lab/Cardiotoxicity/Elisa_24hourCardiotox_20230216.xlsx",
                            sheet = "ExportR", range = "G1:I37",
                            col_names =TRUE)

avgTNNI <- TNNIelisadata %>%
  group_by(indv,Drug) %>%
  summarise(tnni= mean(TNNIavg)) %>%
  mutate(indv=substring(indv,0,2)) %>% #get rid of the -1 on indv to make normalized
  mutate(indv=factor(indv,levels= level_order2)) %>%
  mutate(Drug =case_match(Drug, "Vehicle"~ "Control", .default = Drug))


ggplot(avgTNNI, aes(x=Drug, y=tnni))+
  geom_boxplot(position = "identity")+
  geom_point(aes(col=indv, size=3))+
  geom_signif(comparisons = list(c("Daunorubicin", "Control"),
                                 c("Doxorubicin", "Control"),
                                 c("Epirubicin", "Control"),
                                 c("Mitoxantrone", "Control"),
                                 c("Trastuzumab","Control")),
              test = "t.test",
              map_signif_level = FALSE,step_increase = 0.1,
              textsize = 4)+
  guides(size = "none")+
  ggtitle("Relative troponin I levels released in media")+
  scale_color_brewer(palette = "Dark2",name = "Individual",labels(c(1,2,3,4,5,6)))+
  theme_bw()+
  ylab(" Relative troponin I release")+
  theme(strip.background = element_rect(fill = "transparent")) +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = 15, color = "black"),
    axis.ticks = element_line(size = 1.5),
    axis.text = element_text(size = 9, color = "black", angle = 0),
    strip.text.x = element_text(size = 15, color = "black", face = "bold"))


# LDH 24 hour -------------------------------------------------------------


DA_24_ldh <- matrix(c(1.188,1.222,1.195,1.030,1.074,1.064,1.298,1.282,1.262,
                      1.901,1.975,1.970,3.131,3.246,3.080,1.339,1.438,1.367),
                    ncol =3, nrow =6, byrow =TRUE,
                    dimnames=list(c('87_0.5','79_0.5','75_0.5','77_0.5','78_0.5','71_0.5')))
#ldh24means$Daunorubicin <- rollmean(t(DA_24_ldh),3)
DX_24_ldh <-matrix(c(0.981,0.974,0.978,1.253,1.233,1.292,2.098,2.153,
                     2.114,2.214,2.244,2.239,3.808,3.825,3.735,1.037,1.030,1.030),
                   ncol =3, nrow =6, byrow =TRUE,
                   dimnames=list(c('87_0.5','79_0.5','75_0.5','77_0.5','78_0.5','71_0.5')))
#ldh24means$Doxorubicin <- rollmean(t(DX_24_ldh),3)
EP_24_ldh <-  matrix(c(1.504,1.320,1.469,1.536,1.301,1.531,1.562,1.541,1.558,
                       3.414,3.103,3.236,3.588,3.398,3.611,1.013,0.958,0.991),
                     ncol =3, nrow =6, byrow =TRUE,
                     dimnames=list(c('87_0.5','79_0.5','75_0.5','77_0.5','78_0.5','71_0.5')))
#ldh24means$Epirubicin <- rollmean(t(EP_24_ldh),3)
MT_24_ldh <-  matrix(c(1.508,1.467,1.391,1.493,1.468,1.483,2.010,1.820,1.911,
                       3.089,2.936,2.921,3.623,3.377,3.560,1.222,1.211,1.215),
                     ncol =3, nrow =6, byrow =TRUE,
                     dimnames=list(c('87_0.5','79_0.5','75_0.5','77_0.5','78_0.5','71_0.5')))
#ldh24means$Mitoxantrone <- rollmean(t(MT_24_ldh),3)
TR_24_ldh<-  matrix(c(0.941,0.891,0.953,0.743,0.774,0.812,1.514,1.225,1.252,
                      2.391,1.989,2.172,3.040,2.622,2.613,0.970,0.917,0.895),
                    ncol =3, nrow =6, byrow =TRUE,
                    dimnames=list(c('87_0.5','79_0.5','75_0.5','77_0.5','78_0.5','71_0.5')))
#ldh24means$Trastuzumab <- rollmean(t(TR_24_ldh),3)
VE_24_ldh<-  matrix(c(1.000,1.000,0.977,1.000,1.100,1.096,1.000,0.938,0.951,
                      1.000,1.027,1.038,1.000,1.058,1.062,1.000,1.011,0.975),
                    ncol =3, nrow =6, byrow =TRUE,
                    dimnames=list(c('87_0.5','79_0.5','75_0.5','77_0.5','78_0.5','71_0.5')))
#ldh24means$Vehicle <- rollmean(t(VE_24_ldh),3)
LDH24hstat <- list('VDA'=t.test(VE_24_ldh,DA_24_ldh),
                   'VDX'=t.test(VE_24_ldh,DX_24_ldh),
                   'VEP'=t.test(VE_24_ldh,EP_24_ldh),
                   'VMT'=t.test(VE_24_ldh,MT_24_ldh),
                   'VTR'=t.test(VE_24_ldh,TR_24_ldh),
                   'VVEH'=t.test(VE_24_ldh,VE_24_ldh))
LDH24hstat
t(DA_24_ldh)

mean24ldh <- as.data.frame(rbind(colMeans(t(DA_24_ldh)),
                                 colMeans(t(DX_24_ldh)),
                                 colMeans(t(EP_24_ldh)),
                                 colMeans(t(MT_24_ldh)),
                                 colMeans(t(TR_24_ldh)),
                                 colMeans(t(VE_24_ldh))))
mean24ldh$Drug <- c( "Daunorubicin", "Doxorubicin", "Epirubicin", "Mitoxantrone", "Trastuzumab","Control")  ###add drug name then take out the 0.5 thing
colnames(mean24ldh) <- gsub("_0.5","",colnames(mean24ldh))

mean24ldh <-  mean24ldh %>%
  pivot_longer(.,col=-Drug, names_to = 'indv', values_to = "ldh") %>%
  mutate(indv = factor(indv, levels= level_order2))

ggplot(mean24ldh, aes(x=Drug,y=ldh))+
  geom_boxplot() +
  geom_point(aes(col=indv, size =3))+
  geom_signif(comparisons =list(c("Control","Daunorubicin"),
                                c("Control","Doxorubicin"),
                                c("Control","Epirubicin"),
                                c("Control","Mitoxantrone"),
                                c("Control","Trastuzumab")),
              map_signif_level=FALSE,
              textsize =4,
              tip_length = .1,
              vjust = 0.2, step_increase = 0.1)+
  theme_bw()+
  guides(size = "none")+
  scale_color_brewer(palette = "Dark2")+
  xlab("")+
  ylab("Relative LDH activity ")+
  ggtitle("Relative lactate dehydrogenase release in media")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(0.8)))

# combine data sets tnni and ldh ------------------------------------------
#sets:mean24ldh made above
  #avgTNNI made on the TNNI_ELISA_analysis script

tvl24hour <-  mean24ldh %>%
    full_join(.,avgTNNI, by=c('Drug','indv'))
# LDH 48 hour 0.5 uM ------------------------------------------------------

norm_LDH <- read.csv("data/norm_LDH.csv",row.names=1)##from DRC_analysis.Rmd

norm_LDH$Drug <- factor(norm_LDH$Drug, levels = c("Control", "Daunorubicin", "Doxorubicin", "Epirubicin", "Mitoxantrone", "Trastuzumab"))
norm_LDH %>% filter(Conc==0.5) %>%
  mutate(indv = substring(indv, 0,2)) %>%
  mutate(indv=as.factor(indv)) %>%
  ggplot(., aes(x=Drug,y=norm_val))+
  geom_boxplot() +
  geom_point(aes(col=indv, size =3))+
  geom_signif(comparisons =list(c("Control","Daunorubicin"),
                                c("Control","Doxorubicin"),
                                c("Control","Epirubicin"),
                                c("Control","Mitoxantrone"),
                                c("Control","Trastuzumab")),
              map_signif_level=FALSE,
              textsize =6,
              tip_length = .1,
              vjust = 0.2, step_increase = 0.1)+
  ylab("relative LDH release")+
  theme_bw()+
  ggtitle("48hour lactate dehydrogenase release in media")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(0.8)))


# comparing 48 hour ldh and other types -----------------------------------

norm_24_v_48_ldh <- norm_LDH %>% filter(Conc==0.5) %>%
  mutate(indv=substr(indv,1,2)) %>%
  mutate(indv=as.factor(indv)) %>%
  full_join(., mean24ldh,by=c('Drug','indv'))

ggplot(norm_24_v_48_ldh, aes(x=norm_val, y=ldh))+
  geom_point(aes(col=indv))+
  geom_smooth(method="lm")+
  facet_wrap("Drug")+
  theme_bw()+
  stat_cor(aes(label = after_stat(rr.label)),
           color = "red",
           geom = "label",
           label.y = 4, label.x=1)+
  ylab("48 hour LDH")+
  xlab("24 hour LDH")+
  ggtitle("Relationship between ldh 48 and 24 hours")



ggplot(mean24ldh, aes(x=indv, y=ldh))+geom_point()+ facet_wrap("Drug")
ggplot(avgTNNI), aes(x=indv, y=tnni))+geom_point()+ facet_wrap("Drug")

ggplot(norm_24_v_48_ldh, aes(x=indv, y=norm_val))+geom_point()+ facet_wrap("Drug")


# viability checks --------------------------------------------------------
level_order <- c('71','75','77','78','79','87')
level_order2 <- c('75','87','77','79','78','71')
#write_csv(full_list,"data/Viabilitylistfull.csv")
viafull_list <- read_csv("data/Viabilitylistfull.csv")


via_code <- viafull_list %>%
  #filter(Conc ==0.5) %>%
  mutate(SampleID=substr(SampleID,1,4)) %>%
  mutate(indv= factor(SampleID, levels=c('ind1','ind2', 'ind3', 'ind4', 'ind5','ind6'),
                          labels= level_order2)) %>%
  mutate(Drug=case_match(Drug,"Daun"~"Daunorubicin",
                         "Doxo"~"Doxorubicin",
                         "Epi"~"Epirubicin",
                         "Mito"~"Mitoxantrone",
                         "Tras"~"Trastuzumab",
                         "Veh"~ "Control", .default= Drug))


viability <- via_code %>% dplyr::select(Drug,Name,Conc,Percent,indv) %>%
  group_by(Drug,Conc,indv) %>%
  dplyr::summarize(live=mean(Percent),.groups="drop")


 comparisons_all <- viability %>%full_join(.,norm_24_v_48_ldh, by= c("Drug", "indv")) %>%
 full_join(.,avgTNNI,by= c("Drug", "indv")) %>% select(Drug, indv,live,norm_val,ldh,tnni)
colnames(comparisons_all) <- c('Drug','indv','viable48','ldh_48','ldh_24','tnni_24')

comparisons_allsub <- comparisons_all %>% select('Drug','indv','ldh','tnni')
# 48 hours ldh and viability ----------------------------------------------


ggplot(comparisons_all, aes(x=viable48, y=tnni_24))+
  geom_point(aes(col=indv))+
  geom_smooth(method="lm")+
  facet_wrap("Drug")+
  theme_bw()+
  stat_cor(aes(label = after_stat(rr.label)),
           color = "red",
           geom = "label",
           label.y = 3.75, label.x=0.8)+
  ggtitle("Viability at 48 hours and LDH at 48hours")

# 48 viability and 24 hour ldh --------------------------------------------


ggplot(comparisons_all, aes(x=viable48, y=ldh_24))+
  geom_point(aes(col=indv))+
  geom_smooth(method="lm")+
  facet_wrap("Drug")+
  theme_bw()+
  stat_cor(aes(label = after_stat(rr.label)),
           color = "red",
           geom = "label",
           label.y = 3.75, label.x=.8)+
  ggtitle("Viability at 48 hours and LDH at 24hours")

# 48 viability and 24 tnni ------------------------------------------------

ggplot(comparisons_all, aes(x=tnni_24, y=viable48))+
  geom_point(aes(col=indv))+
  geom_smooth(method="lm")+
  facet_wrap("Drug")+
  theme_bw()+
  stat_cor(aes(label = after_stat(rr.label)),
           color = "red",
           geom = "label",
           label.y = 2.75, label.x=.8)+
  ggtitle("Viability at 48 hours and Troponin I at 24hours")


# normalization using RNA -------------------------------------------------
#using the dataframe tvl24hour made on ldh code block (r combining ldh and tnni data at 24 hours

#write_csv(tvl24hour, "data/tvl24hour.txt")
tvl24hour <- read_csv("data/tvl24hour.txt", show_col_types = FALSE)

tvl24hour <- tvl24hour %>% mutate(indv= factor(indv,levels= level_order2))
#Import the RNA concentratons

all.equal(tvl24hour,comparisons_allsub)
RINsamplelist <- read_csv("data/RINsamplelist.txt")

RNAnormlist <- RINsamplelist %>% mutate(Drug=case_match(Drug,"daunorubicin"~"Daunorubicin",
                                         "doxorubicin"~"Doxorubicin",
                                         "epirubicin"~"Epirubicin",
                                         "mitoxantrone"~"Mitoxantrone",
                                         "trastuzumab"~"Trastuzumab",
                                         "vehicle"~"Control", .default= Drug)) %>%
  filter(time =="24h") %>%
  ungroup() %>%
  dplyr::select(indv,Drug,Conc_ng.ul) %>%
  mutate(indv= factor(indv,levels= level_order2))


ggplot(RNAnormlist, aes(x=indv, y=Conc_ng.ul))+
  geom_point(aes(col=Drug))
ggplot(RNAnormlist, aes(x=Drug, y=Conc_ng.ul))+
  geom_point(aes(col=indv))

RNAnormlist <- RNAnormlist %>%
  full_join(.,tvl24hour,by= c("Drug", "indv")) %>%
  mutate(rldh=ldh/Conc_ng.ul) %>% mutate(rtnni=tnni/Conc_ng.ul)



ggplot(RNAnormlist, aes(x=rldh, y=rtnni))+
  geom_point(aes(col=indv))+
  geom_smooth(method="lm")+
  facet_wrap("Drug", scales="free")+
  theme_bw()+
  scale_color_brewer(palette = "Dark2")+
  #scale_color_manual(values = indv_pal)+
  stat_cor(aes(label = after_stat(rr.label)),
           color = "red",
           geom = "label"
           )+
  ggtitle("Correlation between Troponin I release and LDH activity")



ggplot(RNAnormlist, aes(x=Drug, y=rtnni))+
  geom_boxplot(position = "identity")+
  geom_point(aes(col=indv, size=3))+
  geom_signif(comparisons = list(c("Daunorubicin", "Control"),
                                 c("Doxorubicin", "Control"),
                                 c("Epirubicin", "Control"),
                                 c("Mitoxantrone", "Control"),
                                 c("Trastuzumab","Control")),
              test = "t.test",
              map_signif_level = FALSE,step_increase = 0.1,
              textsize = 4)+
  scale_color_brewer(palette = "Dark2", name = "Individual")+
  ggtitle("Relative troponin I levels released in media after dividing by RNA concentration")+
  theme_bw()+
  guides(size = "none")+
  labs(y = "Relative Troponin I release")+
  #theme(strip.background = element_rect(fill = "transparent")) +
  theme(
    axis.title = element_text(size = 15, color = "black"),
    axis.ticks = element_line(size = 1.5),
    axis.text = element_text(size = 12, color = "black", angle = 0),
    strip.text.x = element_text(size = 15, color = "black", face = "bold"))




ggplot(RNAnormlist, aes(x=Drug, y=ldh))+
  geom_boxplot(position = "identity")+
  geom_point(aes(col=indv, size=3))+
  geom_signif(comparisons = list(c("Daunorubicin", "Control"),
                                 c("Doxorubicin", "Control"),
                                 c("Epirubicin", "Control"),
                                 c("Mitoxantrone", "Control"),
                                 c("Trastuzumab","Control")),
              test = "t.test",
              map_signif_level = FALSE,step_increase = 0.1,
              textsize = 4)+
  scale_color_brewer(palette = "Dark2", name = "Individual")+

  ggtitle("Relative LDH release after dividing by RNA concentration")+
  theme_bw()+
  guides(size = "none")+
  labs(y = "Relative LDH activity", fill = "Individual")+
  #theme(strip.background = element_rect(fill = "transparent")) +
  theme(legend.title = element_text(size = 14),
    axis.title = element_text(size = 15, color = "black"),
    axis.ticks = element_line(size = 1.5),
    axis.text = element_text(size = 12, color = "black", angle = 0),
    strip.text.x = element_text(size = 15, color = "black", face = "bold"))



# Viability graphs --------------------------------------------------------

# viability %>% group_by(Conc) %>%
#   ggplot(., aes(x=as.factor(Conc), y=live, group=as.factor(Conc)))+
#   geom_boxplot()+
#   geom_point(aes(col=indv, size=2))+
#   # #geom_signif(comparisons = list(c("Daunorubicin", "Control"),
#   #                                c("Doxorubicin", "Control"),
#   #                                c("Epirubicin", "Control"),
#   #                                c("Mitoxantrone", "Control"),
#   #             #                    c("Trastuzumab","Control")),
#   #             # test = "t.test",
#               # map_signif_level = FALSE,step_increase = 0.1,
#               # textsize = 6)+
#   scale_color_brewer(palette = "Dark2",name = "Individual", labels = c("6","1","3","5","4","2"))+
#   ggtitle("Relative LDH release")+
#   theme_bw()+
#   facet_wrap(~Drug)
#   guides(size = "none")+
#   labs(y = "Viability Summary", fill = "Individual")+
#   #theme(strip.background = element_rect(fill = "transparent")) +
#   theme(legend.title = element_text(size = 14),
#         axis.title = element_text(size = 15, color = "black"),
#         axis.ticks = element_line(size = 1.5),
#         axis.text = element_text(size = 12, color = "black", angle = 0),
#         strip.text.x = element_text(size = 15, color = "black", face = "bold"))

  LD50long %>%
    ggplot(., mapping = (aes(x = (Treatment), y = log10(mean)))) +
    geom_boxplot()+
    geom_point(aes(
      color = indv,
      size = 5,
      alpha = 0.5)) +
    ggtitle("Observed LD 50s for Cardiomyocytes by Treatment")+
    xlab("Treatment")+
    ylab(bquote('Log'[10]~ 'LD'[50]~'in '*mu~mol))+
    scale_color_brewer(palette = "Dark2",name = "Individual", labels = c("1","2","3","4","5","6"))+
    theme_bw() +
    theme(plot.title = element_text(hjust =0.5, size = 18))+
    guides(alpha ="none", size = "none")
