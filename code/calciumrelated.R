#Calcium analysis
library(tidyverse)
library(ggsignif)
library(readr)

detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base","package:workflowr")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)

}

detachAllPackages()



# load data file ----------------------------------------------------------
clamp_summary <- read.csv("data/Clamp_Summary.csv", row.names=1, col.names =c("Cond","Exp",M) )
level_order2 <- c('75','87','77','79','78','71')
calcium_data <- read_csv("data/DF_Plate_Peak.csv", col_types = cols(...1 = col_skip()))




# Beat rate ---------------------------------------------------------------


calcium_data %>% rename('Treatment'='Condition') %>%
  rename('indv'='Experiment') %>%
  mutate(Drug =case_match(Treatment, "Dau_0.5"~"Daunorubicin",  "Dau_1" ~"Daunorubicin",
                          "Dox_0.5"~"Doxorubicin",  "Dox_1" ~"Doxorubicin",
                          "Epi_0.5"~"Epirubicin",  "Epi_1" ~"Epirubicin",
                          "Mito_0.5"~"Mitoxantrone","Mito_1"~"Mitoxantrone",
                          "Tras_0.5"~"Trastuzumab", "Tras_1"~"Trastuzumab", .default = Treatment)) %>%
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
  geom_boxplot(position = "identity")+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  geom_signif(comparisons = list(c("Daunorubicin", "Control"),
                                 c("Doxorubicin", "Control"),
                                 c("Epirubicin", "Control"),
                                 c("Mitoxantrone", "Control"),
                                 c("Trastuzumab","Control")),
              test = "t.test",
              map_signif_level = TRUE,
              step_increase = 0.1,
              textsize = 4)+
  scale_color_brewer(palette = "Dark2", name = "Individual")+
  ggtitle("Contraction Rate per treatment")+
  theme_bw()+
  guides(size = "none")+
  labs(y = "beat/sec")+
  theme(plot.title = element_text(size=18,hjust = 0.5),
    axis.title = element_text(size = 15, color = "black"),
    axis.ticks = element_line(size = 1.5),
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
                          "Tras_0.5"~"Trastuzumab", "Tras_1"~"Trastuzumab", .default = Treatment)) %>%
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
  geom_signif(comparisons = list(c("Daunorubicin", "Control"),
                                 c("Doxorubicin", "Control"),
                                 c("Epirubicin", "Control"),
                                 c("Mitoxantrone", "Control"),
                                 c("Trastuzumab","Control")),
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



# FWHM --------------------------------------------------------------------



clamp_summary %>%
  rename('Treatment'='Cond') %>%
  rename('indv'='Exp') %>%
  mutate(Drug =case_match(Treatment, "Dau_0.5"~"Daunorubicin",  "Dau_1" ~"Daunorubicin",
                          "Dox_0.5"~"Doxorubicin",  "Dox_1" ~"Doxorubicin",
                          "Epi_0.5"~"Epirubicin",  "Epi_1" ~"Epirubicin",
                          "Mito_0.5"~"Mitoxantrone","Mito_1"~"Mitoxantrone",
                          "Tras_0.5"~"Trastuzumab", "Tras_1"~"Trastuzumab", .default = Treatment)) %>%
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
  geom_signif(comparisons = list(c("Daunorubicin", "Control"),
                                 c("Doxorubicin", "Control"),
                                 c("Epirubicin", "Control"),
                                 c("Mitoxantrone", "Control"),
                                 c("Trastuzumab","Control")),
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
                          "Tras_0.5"~"Trastuzumab", "Tras_1"~"Trastuzumab", .default = Treatment)) %>%
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
  geom_boxplot(position = "identity")+
  geom_point(aes(col=indv, size=2, alpha=0.5))+
  scale_color_brewer(palette = "Dark2", name = "Individual")+
  geom_signif(comparisons = list(c("Daunorubicin", "Control"),
                                 c("Doxorubicin", "Control"),
                                 c("Epirubicin", "Control"),
                                 c("Mitoxantrone", "Control"),
                                 c("Trastuzumab","Control")),
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




facet_wrap(~Drug)
#   guides(size = "none")+
#   labs(y = "Viability Summary", fill = "Individual")+
#   #theme(strip.background = element_rect(fill = "transparent")) +
#   theme(legend.title = element_text(size = 14),
#         axis.title = element_text(size = 15, color = "black"),
#         axis.ticks = element_line(size = 1.5),
#         axis.text = element_text(size = 12, color = "black", angle = 0),
#         strip.text.x = element_text(size = 15, color = "black", face = "bold"))
