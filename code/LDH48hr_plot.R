#adjusting the LDH 48 hour plot
library(readxl)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(zoo)
library(ggplot2)
library(ggsignif)



## Including Plots

# RINsamplelist <- read.csv("~/Ward Lab/Cardiotoxicity/Data/sequencing things/RINsamplelist.csv")
# colnames(RINsamplelist) <- c("Indv","indv","time", "Drug","RIN","Conc_ng.ul")
#  write_csv(RINsamplelist, "data/RINsamplelist.txt")
RINsamplelist <-read_csv("data/RINsamplelist.txt",col_names = TRUE)
factor(RINsamplelist$Drug, levels = c( "daunorubicin", "doxorubicin", "epirubicin", "mitoxantrone", "trastuzumab","control"))
level_order2 <- c('75','87','77','79','78','71')
drug_pal_fact <- c("#8B006D","#DF707E","#F1B72B", "#3386DD","#707031","#41B333")


norm_LDH <- read.csv("data/norm_LDH.csv")
norm_LDH48 <- norm_LDH %>% mutate(Drug=case_match(Drug, "Control" ~"Vehicle",.default= Drug)) %>%
  mutate(Drug=factor(Drug, levels = c( "Daunorubicin", "Doxorubicin", "Epirubicin", "Mitoxantrone", "Trastuzumab","Vehicle"))) %>%
  mutate(indv=substring(indv,0,2)) %>% mutate(indv= factor(indv,levels= level_order2))






# split data out DA----------------------------------------------------------
Daun48 <- norm_LDH48 %>% filter(Drug=="Daunorubicin") %>%
  ggplot(., aes(x = as.factor(Conc), y = norm_val, group = Conc)) +
  geom_boxplot()+#aes(fill = as.factor(Conc))) +

  # scale_fill_manual(values = c("#F3BABD","#EAA2A5","#E28A8E",
  #                              "#DA7277","#CF5157","#C43138","#982124","#6B1210" ),
  #                   aesthetics = c("Conc", "fill")) +
  geom_point(aes(color= indv, alpha= 0.8))+
  scale_color_brewer(palette = "Dark2",name = "Individual")+
  guides(fill = "none", alpha="none", color="none") +
  geom_hline(yintercept = 1,lty =3, lwd=0.8,color="darkgrey") +
  theme_bw() +
  xlab(expression(paste("Drug [", mu, "M]"
  ))) +
  ylab("") +
  labs(fill = "Individual")+
  ylim(0,5)+
  ggtitle("Daunorubicin")+
  theme(plot.title = element_text(hjust = 0.5, size =15, face = "bold"),
        axis.title=element_text(size=10),
        axis.ticks=element_line(linewidth =2),
        axis.text.x = element_text(size = 8, color = "black", face = "bold",angle = 20),
        axis.text.y=element_text(size=8, face = "bold"),
        panel.grid.major = element_line(colour = 'white'),
        panel.border=element_rect(fill = NA, linewidth = 2),
        plot.background = element_rect(fill = "#DDBBCD", colour = NA))


Daun48

# Split Dox ---------------------------------------------------------------

Doxo48 <- norm_LDH48 %>% filter(Drug=="Doxorubicin") %>%
  ggplot(., aes(x = as.factor(Conc), y = norm_val, group = Conc)) +
  geom_boxplot()+
  geom_point(aes(color= indv, alpha= 0.8))+
  scale_color_brewer(palette = "Dark2",name = "Individual")+
  guides(fill = "none", alpha="none", color="none") +
  geom_hline(yintercept = 1,lty =3, lwd=0.8,color="darkgrey") +
  theme_bw() +
  xlab(expression(paste("Drug [", mu, "M]"
  ))) +
 ylab("")+
  labs(fill = "Individual")+
  ylim(0,5)+
  ggtitle("Doxorubicin")+
  theme(plot.title = element_text(hjust = 0.5, size =15, face = "bold"),
        axis.title=element_text(size=10),
        axis.ticks=element_line(linewidth =2),
        axis.text.x = element_text(size = 8, color = "black", face = "bold",angle = 20),
        axis.text.y=element_text(size=8, face = "bold"),
        panel.grid.major = element_line(colour = 'white'),
        panel.border=element_rect(fill = NA, linewidth = 2),
        plot.background = element_rect(fill = "#DDBBCD", colour = NA))
Doxo48
# split Ep ----------------------------------------------------------------


Epi48 <- norm_LDH48 %>% filter(Drug=="Epirubicin") %>%
  ggplot(., aes(x = as.factor(Conc), y = norm_val, group = Conc)) +
  geom_boxplot()+
  geom_point(aes(color= indv, alpha= 0.8))+
  scale_color_brewer(palette = "Dark2",name = "Individual")+
  guides(fill = "none", alpha="none", color="none") +
  geom_hline(yintercept = 1,lty =3, lwd=0.8,color="darkgrey") +
  theme_bw() +
  xlab(expression(paste("Drug [", mu, "M]"
  ))) +
  ylab("") +
  labs(fill = "Individual")+
  ylim(0,5)+
  ggtitle("Epirubicin")+
  theme(plot.title = element_text(hjust = 0.5, size =15, face = "bold"),
        axis.title=element_text(size=10),
        axis.ticks=element_line(linewidth =2),
        axis.text.x = element_text(size = 8, color = "black", face = "bold",angle = 20),
        axis.text.y=element_text(size=8, face = "bold"),
        panel.grid.major = element_line(colour = 'white'),
        panel.border=element_rect(fill = NA, linewidth = 2),
        plot.background = element_rect(fill = "#DDBBCD", colour = NA))

# split mito --------------------------------------------------------------

Mito48 <- norm_LDH48 %>% filter(Drug=="Mitoxantrone") %>%
  ggplot(., aes(x = as.factor(Conc), y = norm_val, group = Conc)) +
  geom_boxplot()+
  geom_point(aes(color= indv, alpha= 0.8))+
  scale_color_brewer(palette = "Dark2",name = "Individual")+
  guides(fill = "none", alpha="none", color="none") +
  geom_hline(yintercept = 1,lty =3, lwd=0.8,color="darkgrey") +
  theme_bw() +
  xlab(expression(paste("Drug [", mu, "M]"
  ))) +
  ylab("") +
  labs(fill = "Individual")+
  ylim(0,5)+
  ggtitle("Mitoxantrone")+
  theme(plot.title = element_text(hjust = 0.5, size =15, face = "bold"),
        axis.title=element_text(size=10),
        axis.ticks=element_line(linewidth =2),
        axis.text.x = element_text(size = 8, color = "black", face = "bold",angle = 20),
        axis.text.y=element_text(size=8, face = "bold"),
        panel.grid.major = element_line(colour = 'white'),
        panel.border=element_rect(fill = NA, linewidth = 2),
        plot.background = element_rect(fill = "#FFD966", colour = NA))
Mito48
# split tras --------------------------------------------------------------



Tras48 <- norm_LDH48 %>% filter(Drug=="Trastuzumab") %>%
  ggplot(., aes(x = as.factor(Conc), y = norm_val, group = Conc)) +
  geom_boxplot()+
  geom_point(aes(color= indv, alpha= 0.8))+
  scale_color_brewer(palette = "Dark2",name = "Individual")+
  guides(fill = "none", alpha="none", color="none") +
  geom_hline(yintercept = 1,lty =3, lwd=0.8,color="darkgrey") +
  theme_bw() +
  xlab(expression(paste("Drug [", mu, "M]"
  ))) +
  ylab("") +
    labs(fill = "Individual")+
  ylim(0,5)+
  ggtitle("Trastuzumab")+
  theme(plot.title = element_text(hjust = 0.5, size =15, face = "bold"),
        axis.title=element_text(size=10),
        axis.ticks=element_line(linewidth =2),
        axis.text.x = element_text(size = 8, color = "black", face = "bold",angle = 20),
        axis.text.y=element_text(size=8, face = "bold"),
        panel.grid.major = element_line(colour = 'white'),
        panel.border=element_rect(fill = NA, linewidth = 2),
        plot.background = element_rect(fill = "#FFD966", colour = NA))




# split Control -----------------------------------------------------------

Con48 <- norm_LDH48 %>% filter(Drug=="Vehicle") %>%
  ggplot(., aes(x = as.factor(Conc), y = norm_val, group = Conc)) +
  geom_boxplot()+
  geom_point(aes(color= indv, alpha= 0.8))+
  scale_color_brewer(palette = "Dark2",name = "Individual")+
  guides(fill = "none", alpha="none", color="none") +
  geom_hline(yintercept = 1,lty =3, lwd=0.8,color="darkgrey") +
  theme_bw() +
  xlab(expression(paste("Drug [", mu, "M]"
  ))) +
  ylab("") +
  labs(fill = "Individual")+
  ylim(0,5)+
  ggtitle("Vehicle")+
  theme(plot.title = element_text(hjust = 0.5, size =15, face = "bold"),
        axis.title=element_text(size=10),
        axis.ticks=element_line(linewidth =2),
        axis.text.x = element_text(size = 8, color = "black", face = "bold",angle = 20),
        axis.text.y=element_text(size=8, face = "bold"),
        panel.grid.major = element_line(colour = 'white'),
        panel.border=element_rect(fill = NA, linewidth = 2),
        plot.background = element_rect(fill = "transparent", colour = NA))

Con48
# all together ------------------------------------------------------------





library(cowplot)
## take a legend from a plot
#legend_b <- ggpubr::get_legend(leg)
#
# ggpubr::as_ggplot(legend_b)
# legend_b

# plan1 <- plot_grid(daplot,dxplot,epplot,NULL,mtplot, trplot,veplot,legend_b,ncol =4)
plan48ldh <-  plot_grid(Daun48,Doxo48,
                        Epi48, Mito48, Tras48, Con48,ncol =3)

#drcfinal<- plot_grid(plan1, ncol =1, rel_heights = c(1,0.1))
save_plot('output/plan48ldh.png', plan48ldh)
#save_plot('individual-legenddark2.png',ggpubr::as_ggplot(legend_b))
plan48ldh
