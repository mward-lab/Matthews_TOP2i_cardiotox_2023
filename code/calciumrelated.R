#Calcium analysis
library(tidyverse)
library(ggsignif)
library(readr)
library(ggalt)

detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base","package:workflowr")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)

}

detachAllPackages()



# load data file ----------------------------------------------------------
clamp_summary <- read.csv("data/Clamp_Summary.csv", row.names=1)

level_order2 <- c('75','87','77','79','78','71')
calcium_data <- read_csv("data/DF_Plate_Peak.csv", col_types = cols(...1 = col_skip()))

drug_palexpand <- c("#41B333","#8B006D","#DF707E","#F1B72B", "#3386DD","#707031","purple3","darkgreen", "darkblue")
#named colors: dark pink,Red,yellow,blue, dark grey, green(green is always control, may need to move pal around)


# Beat rate ---------------------------------------------------------------


calcium_data %>% rename('Treatment'='Condition') %>%
  rename('indv'='Experiment') %>%
  mutate(Drug =case_match(Treatment, "Dau_0.5"~"Daunorubicin",  "Dau_1" ~"Daunorubicin",
                          "Dox_0.5"~"Doxorubicin",  "Dox_1" ~"Doxorubicin",
                          "Epi_0.5"~"Epirubicin",  "Epi_1" ~"Epirubicin",
                          "Mito_0.5"~"Mitoxantrone","Mito_1"~"Mitoxantrone",
                          "Tras_0.5"~"Trastuzumab", "Tras_1"~"Trastuzumab","Control" ~"Vehicle", .default = Treatment)) %>%
  mutate(Drug=factor(Drug, levels = c("Vehicle", "Daunorubicin", "Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab"))) %>%
  mutate(Conc =case_match(Treatment, "Dau_0.5"~"0.5",  "Dau_1" ~"1.0",
                          "Dox_0.5"~"0.5",  "Dox_1" ~"1.0",
                          "Epi_0.5"~"0.5",  "Epi_1" ~"1.0",
                          "Mito_0.5"~"0.5","Mito_1"~"1.0",
                          "Tras_0.5"~"0.5", "Tras_1"~"1.0",'Control'~'0', .default = Treatment)) %>%
  dplyr::select(Drug,Conc,indv,Rate) %>%#Peak_variance,Ave_FF0,
  mutate(indv=substr(indv,1,2)) %>%
  mutate(indv=factor(indv, levels = level_order2)) %>%

  mutate(contrl= 0.383) %>%
  mutate(norm_rate=Rate/contrl) %>%
  filter(Conc==0| Conc==0.5) %>%
  ggplot(., aes(x=Drug, y=Rate))+
  geom_boxplot(position = "identity", fill= drug_pal)+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  guides(alpha= "none")+
  geom_signif(comparisons = list(c("Daunorubicin", "Vehicle"),
                                 c("Doxorubicin", "Vehicle"),
                                 c("Epirubicin", "Vehicle"),
                                 c("Mitoxantrone", "Vehicle"),
                                 c("Trastuzumab","Vehicle")),
              test = "t.test",
              map_signif_level = TRUE,
              step_increase = 0.1,
              textsize = 4)+
  scale_color_brewer(palette = "Dark2", name = "Individual")+
  ggtitle("Contraction rate")+
  theme_classic()+
  guides(size = "none",colour = guide_legend(override.aes = list(size=4, alpha= 0.5)))+
  labs(y = "avg. beats/sec")+
  theme(plot.title = element_text(size=18,hjust = 0.5),
    axis.title = element_text(size = 15, color = "black"),
    axis.ticks = element_line(linewidth = 1.5),
    axis.line = element_line(linewidth = 1.5),
    axis.text = element_text(size = 12, color = "black", angle = 0),
    strip.text.x = element_text(size = 15, color = "black", face = "bold"))


# Amplitude ---------------------------------------------------------------


clamp_summary %>%
  rename('Treatment'='Cond') %>%
  rename('indv'='Exp') %>%
  mutate(Drug =case_match(Treatment, "Dau_0.5"~"Daunorubicin",  "Dau_1" ~"Daunorubicin",
                          "Dox_0.5"~"Doxorubicin",  "Dox_1" ~"Doxorubicin",
                          "Epi_0.5"~"Epirubicin",  "Epi_1" ~"Epirubicin",
                          "Mito_0.5"~"Mitoxantrone","Mito_1"~"Mitoxantrone",
                          "Tras_0.5"~"Trastuzumab", "Tras_1"~"Trastuzumab","Control" ~"Vehicle" ,.default = Treatment)) %>%
  mutate(Drug=factor(Drug, levels = c("Vehicle", "Daunorubicin", "Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab"))) %>%
  mutate(Conc =case_match(Treatment, "Dau_0.5"~"0.5",  "Dau_1" ~"1.0",
                          "Dox_0.5"~"0.5",  "Dox_1" ~"1.0",
                          "Epi_0.5"~"0.5",  "Epi_1" ~"1.0",
                          "Mito_0.5"~"0.5","Mito_1"~"1.0",
                          "Tras_0.5"~"0.5", "Tras_1"~"1.0",'Control'~'0', .default = Treatment)) %>%
  rename(c('Mean_Amplitude'='R1S1_Mean..a.u..',
           'Rise_Slope'='R1S1_Rise_Slope..a.u..ms.',
           'FWHM'='R1S1_Half_Width..ms.',
           'Decay_Slope'='R1S1_Decay_Slope..a.u..ms.',
           'Decay_Time'='Decay_Time..ms.',
           'Rise_Time'='R1S1_Rise_Time')) %>%
  mutate(indv=substr(indv,1,2)) %>%
  mutate(indv=factor(indv, levels = level_order2)) %>%
  filter(Conc==0| Conc==0.5) %>%
  dplyr::select(Drug,Conc,indv,Mean_Amplitude,FWHM,Rise_Slope,Decay_Slope) %>%
  ggplot(.,aes(Drug,Mean_Amplitude))+
  geom_boxplot(position = "identity", fill= drug_pal)+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  guides(alpha= "none",size="none")+
  scale_color_brewer(palette = "Dark2", name = "Individual",label=c("2","6","5"))+
  geom_signif(comparisons = list(c("Daunorubicin", "Vehicle"),
                                 c("Doxorubicin", "Vehicle"),
                                 c("Epirubicin", "Vehicle"),
                                 c("Mitoxantrone", "Vehicle"),
                                 c("Trastuzumab","Vehicle")),
              test = "t.test",
              map_signif_level = TRUE,
              step_increase = 0.1,
              textsize = 4)+
  labs(y = "a.u.")+
  # ggtitle(expression(paste("Mean amplitude at 0.5 ",mu,"M")))+
  ggtitle("Mean amplitude")+
  theme_classic()+
  theme(plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text = element_text(size = 12, color = "black", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))



# FWHM --------------------------------------------------------------------



clamp_summary %>%
  rename('Treatment'='Cond') %>%
  rename('indv'='Exp') %>%
  mutate(Drug =case_match(Treatment, "Dau_0.5"~"Daunorubicin",  "Dau_1" ~"Daunorubicin",
                          "Dox_0.5"~"Doxorubicin",  "Dox_1" ~"Doxorubicin",
                          "Epi_0.5"~"Epirubicin",  "Epi_1" ~"Epirubicin",
                          "Mito_0.5"~"Mitoxantrone","Mito_1"~"Mitoxantrone",
                          "Tras_0.5"~"Trastuzumab", "Tras_1"~"Trastuzumab","Control" ~"Vehicle" ,.default = Treatment)) %>%
  mutate(Drug=factor(Drug, levels = c("Vehicle", "Daunorubicin", "Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab"))) %>%
  mutate(Conc =case_match(Treatment, "Dau_0.5"~"0.5",  "Dau_1" ~"1.0",
                          "Dox_0.5"~"0.5",  "Dox_1" ~"1.0",
                          "Epi_0.5"~"0.5",  "Epi_1" ~"1.0",
                          "Mito_0.5"~"0.5","Mito_1"~"1.0",
                          "Tras_0.5"~"0.5", "Tras_1"~"1.0",'Control'~'0', .default = Treatment)) %>%
  rename(c('Mean_Amplitude'='R1S1_Mean..a.u..',
           'Rise_Slope'='R1S1_Rise_Slope..a.u..ms.',
           'FWHM'='R1S1_Half_Width..ms.',
           'Decay_Slope'='R1S1_Decay_Slope..a.u..ms.',
           'Decay_Time'='Decay_Time..ms.',
           'Rise_Time'='R1S1_Rise_Time')) %>%
  mutate(indv=substr(indv,1,2)) %>%
  mutate(indv=factor(indv, levels = level_order2)) %>%
  filter(Conc==0| Conc==0.5) %>%
  dplyr::select(Drug,Conc,indv,Mean_Amplitude,FWHM,Rise_Slope,Decay_Slope) %>%
  ggplot(.,aes(Drug,Mean_Amplitude))+
  geom_boxplot(position = "identity")+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  scale_color_brewer(palette = "Dark2", name = "Individual")+
  geom_signif(comparisons = list(c("Daunorubicin", "Vehicle"),
                                 c("Doxorubicin", "Vehicle"),
                                 c("Epirubicin", "Vehicle"),
                                 c("Mitoxantrone", "Vehicle"),
                                 c("Trastuzumab","Vehicle")),
              test = "t.test",
              map_signif_level = TRUE,
              step_increase = 0.1,
              textsize = 4)+
  ggtitle(expression(paste("Mean Amplitude at 0.5 ",mu,"M")))+
  theme_bw()+
  guides(size = "none")+
  labs(y = "a.u.")+
  theme(plot.title=element_text(size= 18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(size = 1.5),
        axis.text = element_text(size = 12, color = "black", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))







# rising slope ------------------------------------------------------------


clamp_summary %>%
  rename('Treatment'='Cond') %>%
  rename('indv'='Exp') %>%
  mutate(Drug =case_match(Treatment, "Dau_0.5"~"Daunorubicin",  "Dau_1" ~"Daunorubicin",
                          "Dox_0.5"~"Doxorubicin",  "Dox_1" ~"Doxorubicin",
                          "Epi_0.5"~"Epirubicin",  "Epi_1" ~"Epirubicin",
                          "Mito_0.5"~"Mitoxantrone","Mito_1"~"Mitoxantrone",
                          "Tras_0.5"~"Trastuzumab", "Tras_1"~"Trastuzumab","Control" ~"Vehicle" ,.default = Treatment)) %>%
  mutate(Drug=factor(Drug, levels = c("Vehicle", "Daunorubicin", "Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab"))) %>%
  mutate(Conc =case_match(Treatment, "Dau_0.5"~"0.5",  "Dau_1" ~"1.0",
                          "Dox_0.5"~"0.5",  "Dox_1" ~"1.0",
                          "Epi_0.5"~"0.5",  "Epi_1" ~"1.0",
                          "Mito_0.5"~"0.5","Mito_1"~"1.0",
                          "Tras_0.5"~"0.5", "Tras_1"~"1.0",'Control'~'0', .default = Treatment)) %>%
  rename(c('Mean_Amplitude'='R1S1_Mean..a.u..',
           'Rise_Slope'='R1S1_Rise_Slope..a.u..ms.',
           'FWHM'='R1S1_Half_Width..ms.',
           'Decay_Slope'='R1S1_Decay_Slope..a.u..ms.',
           'Decay_Time'='Decay_Time..ms.',
           'Rise_Time'='R1S1_Rise_Time')) %>%
  mutate(indv=substr(indv,1,2)) %>%
  mutate(indv=factor(indv, levels = level_order2)) %>%
  filter(Conc==0| Conc==0.5) %>%
  dplyr::select(Drug,Conc,indv,Mean_Amplitude,FWHM,Rise_Slope,Decay_Slope) %>%
  ggplot(.,aes(Drug,Rise_Slope))+
  geom_boxplot(position = "identity", fill= drug_pal)+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  guides(alpha="none", indv = "none", size = "none")+
  scale_color_brewer(palette = "Dark2", name = "Individual")+
  geom_signif(comparisons = list(c("Daunorubicin", "Vehicle"),
                                 c("Doxorubicin", "Vehicle"),
                                 c("Epirubicin", "Vehicle"),
                                 c("Mitoxantrone", "Vehicle"),
                                 c("Trastuzumab","Vehicle")),
              test = "t.test",
              map_signif_level = TRUE,
              step_increase = 0.1,
              textsize = 4)+
  ggtitle(expression(paste("Rising slope")))+
  labs(y = "a.u.")+
  theme_classic()+
  theme(plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text = element_text(size = 12, color = "black", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))




# Decay slope -------------------------------------------------------------
clamp_summary %>%
  rename('Treatment'='Cond') %>%
  rename('indv'='Exp') %>%
  mutate(Drug =case_match(Treatment, "Dau_0.5"~"Daunorubicin",  "Dau_1" ~"Daunorubicin",
                          "Dox_0.5"~"Doxorubicin",  "Dox_1" ~"Doxorubicin",
                          "Epi_0.5"~"Epirubicin",  "Epi_1" ~"Epirubicin",
                          "Mito_0.5"~"Mitoxantrone","Mito_1"~"Mitoxantrone",
                          "Tras_0.5"~"Trastuzumab", "Tras_1"~"Trastuzumab","Control" ~"Vehicle" ,.default = Treatment)) %>%
  mutate(Drug=factor(Drug, levels = c("Vehicle", "Daunorubicin", "Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab"))) %>%
  mutate(Conc =case_match(Treatment, "Dau_0.5"~"0.5",  "Dau_1" ~"1.0",
                          "Dox_0.5"~"0.5",  "Dox_1" ~"1.0",
                          "Epi_0.5"~"0.5",  "Epi_1" ~"1.0",
                          "Mito_0.5"~"0.5","Mito_1"~"1.0",
                          "Tras_0.5"~"0.5", "Tras_1"~"1.0",'Control'~'0', .default = Treatment)) %>%
  rename(c('Mean_Amplitude'='R1S1_Mean..a.u..',
           'Rise_Slope'='R1S1_Rise_Slope..a.u..ms.',
           'FWHM'='R1S1_Half_Width..ms.',
           'Decay_Slope'='R1S1_Decay_Slope..a.u..ms.',
           'Decay_Time'='Decay_Time..ms.',
           'Rise_Time'='R1S1_Rise_Time')) %>%
  mutate(indv=substr(indv,1,2)) %>%
  mutate(indv=factor(indv, levels = level_order2)) %>%
  filter(Conc==0| Conc==0.5) %>%
  dplyr::select(Drug,Conc,indv,Mean_Amplitude,FWHM,Rise_Slope,Decay_Slope) %>%
  ggplot(.,aes(Drug,Decay_Slope))+
  geom_boxplot(position = "identity", fill= drug_pal)+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  guides(alpha="none", indv = "none", size = "none")+
  scale_color_brewer(palette = "Dark2", name = "Individual")+
  geom_signif(comparisons = list(c("Daunorubicin", "Vehicle"),
                                 c("Doxorubicin", "Vehicle"),
                                 c("Epirubicin", "Vehicle"),
                                 c("Mitoxantrone", "Vehicle"),
                                 c("Trastuzumab","Vehicle")),
              test = "t.test",
              map_signif_level = TRUE,
              step_increase = 0.1,
              textsize = 4)+
  ggtitle(expression(paste("Decay slope ")))+
  labs(y = "a.u.")+
  theme_classic()+
  theme(plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text = element_text(size = 12, color = "black", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))
# playing with rising/decay -----------------------------------------------

clamp_summary %>%
  rename('Treatment'='Cond') %>%
  rename('indv'='Exp') %>%
  mutate(Drug =case_match(Treatment, "Dau_0.5"~"Daunorubicin",  "Dau_1" ~"Daunorubicin",
                          "Dox_0.5"~"Doxorubicin",  "Dox_1" ~"Doxorubicin",
                          "Epi_0.5"~"Epirubicin",  "Epi_1" ~"Epirubicin",
                          "Mito_0.5"~"Mitoxantrone","Mito_1"~"Mitoxantrone",
                          "Tras_0.5"~"Trastuzumab", "Tras_1"~"Trastuzumab","Control" ~"Vehicle" .default = Treatment)) %>%
  mutate(Conc =case_match(Treatment, "Dau_0.5"~"0.5",  "Dau_1" ~"1.0",
                          "Dox_0.5"~"0.5",  "Dox_1" ~"1.0",
                          "Epi_0.5"~"0.5",  "Epi_1" ~"1.0",
                          "Mito_0.5"~"0.5","Mito_1"~"1.0",
                          "Tras_0.5"~"0.5", "Tras_1"~"1.0",'Vehicle'~'0', .default = Treatment)) %>%
  rename(c('Mean_Amplitude'='R1S1_Mean..a.u..',
           'Rise_Slope'='R1S1_Rise_Slope..a.u..ms.',
           'FWHM'='R1S1_Half_Width..ms.',
           'Decay_Slope'='R1S1_Decay_Slope..a.u..ms.',
           'Decay_Time'='Decay_Time..ms.',
           'Rise_Time'='R1S1_Rise_Time')) %>%
  mutate(indv=substr(indv,1,2)) %>%
  mutate(indv=factor(indv, levels = level_order2)) %>%
  filter(Conc==0| Conc==0.5) %>%
  dplyr::select(Drug,Conc,indv,Mean_Amplitude,FWHM,Rise_Slope,Decay_Slope) %>%
  mutate(fraction= Rise_Slope/-Decay_Slope) %>%
  ggplot(.,aes(Drug,fraction))+
  geom_boxplot(position = "identity")+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  scale_color_brewer(palette = "Dark2", name = "Individual")+
  geom_signif(comparisons = list(c("Daunorubicin", "Vehicle"),
                                 c("Doxorubicin", "Vehicle"),
                                 c("Epirubicin", "Vehicle"),
                                 c("Mitoxantrone", "Vehicle"),
                                 c("Trastuzumab","Vehicle")),
              test = "t.test",
              map_signif_level = TRUE,
              step_increase = 0.1,
              textsize = 4)+
  ggtitle(expression(paste("Ratio of CA2+ exchange Rise/Decay at 0.5 ",mu,"M")))+
  theme_bw()+
  guides(size = "none")+
  labs(y = "a.u.")+
  theme(plot.title=element_text(size= 18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(size = 1.5),
        axis.text = element_text(size = 12, color = "black", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))

# k means -----------------------------------------------------------------

k_means <- read.csv("data/K_cluster_kisthree.csv")

k_means %>% mutate(Drug =case_match(Drug_Name, "Dau_0.5"~"Daunorubicin",  "Dau_0.5.1" ~"Daunorubicin","Dau_0.5.2" ~"Daunorubicin",
                                    "Dox_0.5"~"Doxorubicin",  "Dox_0.5.1" ~"Doxorubicin","Dox_0.5.2" ~"Doxorubicin",
                                    "Epi_0.5"~"Epirubicin",  "Epi_0.5.1"~"Epirubicin","Epi_0.5.2"~"Epirubicin",
                                    "Mito_0.5"~"Mitoxantrone","Mito_0.5.1"~"Mitoxantrone","Mito_0.5.2"~"Mitoxantrone",
                                    "Tras_0.5"~"Trastuzumab", "Tras_0.5.1"~"Trastuzumab", "Tras_0.5.2"~"Trastuzumab",
                                    "Control.1"~"Vehicle","Control.2"~"Vehicle","Control"~"Vehicle",.default = Drug_Name)) %>%
  ggplot(., aes(x=PC1, y=PC2, col=factor(Cluster), shape= Cluster))+
  geom_point(size=5)+
  scale_color_manual(values=c("purple3","darkgreen", "darkblue"))


k_means %>% mutate(Drug =case_match(Drug_Name, "Dau_0.5"~"Daunorubicin",  "Dau_0.5.1" ~"Daunorubicin","Dau_0.5.2" ~"Daunorubicin",
                                    "Dox_0.5"~"Doxorubicin",  "Dox_0.5.1" ~"Doxorubicin","Dox_0.5.2" ~"Doxorubicin",
                                    "Epi_0.5"~"Epirubicin",  "Epi_0.5.1"~"Epirubicin","Epi_0.5.2"~"Epirubicin",
                                    "Mito_0.5"~"Mitoxantrone","Mito_0.5.1"~"Mitoxantrone","Mito_0.5.2"~"Mitoxantrone",
                                    "Tras_0.5"~"Trastuzumab", "Tras_0.5.1"~"Trastuzumab", "Tras_0.5.2"~"Trastuzumab",
                                    "Control.1"~"Vehicle","Control.2"~"Vehicle","Control"~"Vehicle",.default = Drug_Name)) %>%
  mutate(Drug=factor(Drug, levels = c("Vehicle", "Daunorubicin", "Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab"))) %>%
  mutate(Cluster=factor(Cluster)) %>%
  ggplot(., aes(x=PC1, y=PC2, col= Drug))+
  geom_point(shape = 17,size=7)+
  #stat_ellipse +f
 scale_color_manual(values=drug_pal)+
  #scale_shape_discrete(name = "Cluster")+
 # geom_encircle(data= Cluster,s_shape=0.5, expand =0)+
  geom_encircle(aes(group=Cluster))+
  annotate("text", label = c("Cluster 1","Cluster 2", "Cluster 3"), x = c(-2,0,1.5),y=c(-0.5,0,0.5))+
  ggtitle(expression("PCA of Ca"^("2+")~"data with K-means clustering"))+
  theme_bw()+
  #guides(size = "none")+
  labs(x = "PC 1 (54 %)",y = "PC 2 (34%)")+
  theme(plot.title=element_text(size= 14,hjust = 0.5),
        axis.title = element_text(size = 10, color = "black"),
        axis.ticks = element_line(size = 1.5),
        axis.text = element_text(size = 12, color = "black", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))
