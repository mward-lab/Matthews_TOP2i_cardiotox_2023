
##r script to compare gene lists in common with Top2Bi genes and Top2Bi from edgeR data
# load libraries ---------------------------------------------------------

library(tidyverse)
library(limma)
library(readr)
library(BiocGenerics)
library(gridExtra)
library(VennDiagram)
library(kableExtra)
library(scales)
library(ggVennDiagram)
library(cowplot)
library(RColorBrewer)
library(gprofiler2)

# load data sets ----------------------------------------------------------

##data set of cormotif Top2Bi expressed genes
DEG_cormotif <- readRDS("data/DEG_cormotif.RDS")
motif1_NR <- DEG_cormotif$motif1_NR
motif3_TI <- DEG_cormotif$motif3_TI
motif4_LR <- DEG_cormotif$motif4_LR
motif5_ER <- DEG_cormotif$motif5_ER

backGL <- read.csv("data/backGL.txt")
NRresp <- read_csv("data/cormotif_NRset.txt")


##data set of AConly response genes from edgeR analysis (aka 4400)
file.names <- list.files(path = "data/", pattern = "sig*", ignore.case = TRUE,full.names = TRUE)

filenameonly <- read_csv("data/filenameonly.txt")
#loop through the list of files and make a separate dataframe for each file under the 'real'  name of the data set
for (k in 1:length(file.names)){

  assign(paste0(filenameonly$x[k]) , read.csv(file.names[k]))
}

##rename the columns to the previous names

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

# venn diagram overlap retrieval ------------------------------------------


total24 <-list(sigVDA24$ENTREZID,sigVDX24$ENTREZID,sigVEP24$ENTREZID,sigVMT24$ENTREZID)

list24totvenn <- get.venn.partitions(total24)
DDEresp<- list24totvenn$..values..[[9]]
#write.csv(DDEresp,"data/DDEresp_list.csv")
DDEMresp<- list24totvenn$..values..[[1]]

# venn diagrams of data sets -----------------------------------------------

# ERmotif and AC resp -----------------------------------------------------


list_methodcomp <- list(as.numeric(motif5_ER),as.numeric(DDEresp))

list24tovennERAC <- get.venn.partitions(list_methodcomp)


ggVennDiagram(list_methodcomp,
              category.names = c("ERmotif","AC response"),
              show_intersect = FALSE,
              set_color = "black",
              catagory_size = c(6,6,6,6),
              label = "both",
              label_percent_digit = 1,
              label_size = 3,
              label_alpha = 0,
              label_color = "black",
              edge_lty = "solid", set_size = )+
  scale_x_continuous(expand = expansion(mult = .3))+
  scale_y_continuous(expand = expansion(mult = .2))+
  #scale_color_brewer(palette = "Dark2",name = "Individual",)
  scale_fill_distiller(palette="RdYlBu")+
  labs(title = "ER and AC overlaps", caption = "n =genes")+
  theme(plot.title = element_text(size = rel(1.6), hjust = 0.5, vjust =1))

# LR overlap --------------------------------------------------------------

list_methodcomp1 <- list(as.numeric(motif4_LR),as.numeric(DDEresp))
list24totvennLRAC <- get.venn.partitions(list_methodcomp1)
ggVennDiagram(list_methodcomp1,
              category.names = c("LRmotif","AC response"),
              show_intersect = FALSE,
              set_color = "black",
              catagory_size = c(6,6,6,6),
              label = "both",
              label_percent_digit = 1,
              label_size = 3,
              label_alpha = 0,
              label_color = "black",
              edge_lty = "solid", set_size = )+
  scale_x_continuous(expand = expansion(mult = .3))+
  scale_y_continuous(expand = expansion(mult = .2))+
  #scale_color_brewer(palette = "Dark2",name = "Individual",)
  scale_fill_distiller(palette="RdYlBu")+
  labs(title = "LR motif and AC degs", caption = "n = 8887 genes")+
  theme(plot.title = element_text(size = rel(1.6), hjust = 0.5, vjust =1))

# TI motif and ac ---------------------------------------------------------

list_methodcompti <- list(as.numeric(motif3_TI),as.numeric(DDEresp))
list24totvennTIAC <- get.venn.partitions(list_methodcompti)

ggVennDiagram(list_methodcompti,
              category.names = c("TI motif","AC response"),
              show_intersect = FALSE,
              set_color = "black",
              catagory_size = c(6,6,6,6),
              label = "both",
              label_percent_digit = 1,
              label_size = 3,
              label_alpha = 0,
              label_color = "black",
              edge_lty = "solid", set_size = )+
  scale_x_continuous(expand = expansion(mult = .3))+
  scale_y_continuous(expand = expansion(mult = .2))+
  #scale_color_brewer(palette = "Dark2",name = "Individual",)
  scale_fill_distiller(palette="RdYlBu")+
  labs(title = "TI motif and AC degs", caption = "n = 8887 genes")+
  theme(plot.title = element_text(size = rel(1.6), hjust = 0.5, vjust =1))





# top2bi comparison both sets ----------------------------------------------



list_methodcomptop1 <- list(as.numeric(motif5_ER),as.numeric(DDEresp))
list24totvennERTop <- get.venn.partitions(list_methodcomptop1)


ggVennDiagram(list_methodcomptop1,
              category.names = c("ERmotif","Top2bi"),
              show_intersect = FALSE,
              set_color = "black",
              catagory_size = c(6,6,6,6),
              label = "both",
              label_percent_digit = 1,
              label_size = 3,
              label_alpha = 0,
              label_color = "black",
              edge_lty = "solid", set_size = )+
  scale_x_continuous(expand = expansion(mult = .3))+
  scale_y_continuous(expand = expansion(mult = .2))+
  #scale_color_brewer(palette = "Dark2",name = "Individual",)
  scale_fill_distiller(palette="RdYlBu")+
  labs(title = "ER and top2bi overlaps", caption = "n =genes")+
  theme(plot.title = element_text(size = rel(1.6), hjust = 0.5, vjust =1))

# LR overlap DDEM --------------------------------------------------------------
list24totvennLRTop <- get.venn.partitions(list_methodcomptop2)

list_methodcomptop2 <- list(as.numeric(motif4_LR),as.numeric(DDEMresp))

ggVennDiagram(list_methodcomptop2,
              category.names = c("LRmotif","top2bi response"),
              show_intersect = FALSE,
              set_color = "black",
              catagory_size = c(6,6,6,6),
              label = "both",
              label_percent_digit = 1,
              label_size = 3,
              label_alpha = 0,
              label_color = "black",
              edge_lty = "solid", set_size = )+
  scale_x_continuous(expand = expansion(mult = .3))+
  scale_y_continuous(expand = expansion(mult = .2))+
  #scale_color_brewer(palette = "Dark2",name = "Individual",)
  scale_fill_distiller(palette="RdYlBu")+
  labs(title = "LR motif and top2bi degs", caption = "n = 8887 genes")+
  theme(plot.title = element_text(size = rel(1.6), hjust = 0.5, vjust =1))

# TI motif and ac ---------------------------------------------------------
list24totvennTITop <- get.venn.partitions(list_methodcomptop3)
list_methodcomptop3 <- list(as.numeric(motif3_TI),as.numeric(DDEMresp))

ggVennDiagram(list_methodcomptop3,
              category.names = c("TI motif","top2bi response"),
              show_intersect = FALSE,
              set_color = "black",
              catagory_size = c(6,6,6,6),
              label = "both",
              label_percent_digit = 1,
              label_size = 3,
              label_alpha = 0,
              label_color = "black",
              edge_lty = "solid", set_size = )+
  scale_x_continuous(expand = expansion(mult = .3))+
  scale_y_continuous(expand = expansion(mult = .2))+
  #scale_color_brewer(palette = "Dark2",name = "Individual",)
  scale_fill_distiller(palette="RdYlBu")+
  labs(title = "TI motif and Top2Bi degs", caption = "n = genes")+
  theme(plot.title = element_text(size = rel(1.6), hjust = 0.5, vjust =1))


# GO analysis of AC and Top2Bi Overlaps -----------------------------------

# AC only 24 hours vs top2bi LR genes -------------------------------------
list24tovennERAC
list24totvennLRAC
list24totvennTIAC
early_AConly <- list24tovennERAC$..values..[[1]]
late_AConly <- list24totvennLRAC$..values..[[1]]
TI_AConly <- list24totvennTIAC$..values..[[1]]


LR24ACgost <- gost(query = late_AConly,
                      organism = "hsapiens",
                      ordered_query = TRUE,
                      domain_scope = "custom",
                      measure_underrepresentation = FALSE,
                      evcodes = FALSE,
                      user_threshold = 0.05,
                      correction_method = c("fdr"),
                      custom_bg = backGL$ENTREZID,
                      sources=c("GO:BP","KEGG"))

LR24ACtable <- LR24ACgost$result %>%
  dplyr::select(c(source, term_id, term_name,intersection_size, term_size, p_value))
##only one kegg term no GObp out of 510 overlap


# LR response and TOP2bi specific sets ------------------------------------

LR24topval <- list24totvennLRTop$..values..[[1]]
LR24topgost <- gost(query = LR24topval,
                      organism = "hsapiens",
                      ordered_query = TRUE,
                      domain_scope = "custom",
                      measure_underrepresentation = FALSE,
                      evcodes = FALSE,
                      user_threshold = 0.05,
                      correction_method = c("fdr"),
                      custom_bg = backGL$ENTREZID,
                      sources=c("GO:BP","KEGG"))

LR24toptable <- LR24topgost$result %>%
  dplyr::select(c(source, term_id, term_name,intersection_size, term_size, p_value))



LR24toptable%>%
  dplyr::filter(source=="GO:BP") %>%
  dplyr::select(p_value,term_name,intersection_size) %>%
  slice_min(., n=10 ,order_by=p_value) %>%
  mutate(log_val = -log10(p_value)) %>%
  # slice_max(., n=10,order_by = p_value) %>%
  ggplot(., aes(x = log_val, y =reorder(term_name,p_value), col= intersection_size)) +
  geom_point(aes(size = intersection_size)) +
  scale_y_discrete(labels = wrap_format(30))+
  guides(col="none", size=guide_legend(title = "# of intersected \n terms"))+
  ggtitle("Mitoxantrone specific 24 hour top GO:BP terms") +
  xlab(expression(" -"~log[10]~"(adj p-value)"))+
  ylab("GO: BP term")+
  theme_bw()+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text = element_text(size = 10, color = "black", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))



# knowles interface -------------------------------------------------------

##here is where I compare eqtl enrichment from knowles data with the 4400 and 952 sets (dde/ddem)
total24 <-list(sigVDA24$ENTREZID,sigVDX24$ENTREZID,sigVEP24$ENTREZID,sigVMT24$ENTREZID)

list24totvenn <- get.venn.partitions(total24)
DDEresp<- list24totvenn$..values..[[9]]
#write.csv(DDEresp,"data/DDEresp_list.csv")
DDEMresp<- list24totvenn$..values..[[1]]

# now we see how these sets over lap with the K4 and K5 sets

##24hours no response
##No response gene set
NoResp$ENTREZID
DDEresp_genes <- toplist24hours %>%
  dplyr::filter(ENTREZID%in%DDEresp)

Overlapk4DDE <- DDEresp_genes %>%
  dplyr::filter(ENTREZID%in%knowles4) %>%
  group_by(id) %>%

 count()
intersect(NoResp$ENTREZID,knowles4) ##find noresp for k4
NR4<- tibble_row( id = "No Response", n = 248)
Overlapk4 <- merge(Overlapk4,NR4, all =TRUE)[,1:2]




Overlapk5 <- toplist24hours %>%
  dplyr::filter(ENTREZID%in%knowles5) %>%
  group_by(id) %>%
  filter(adj.P.Val<0.1) %>%
  count()

NRk5 <- tibble_row( id = "No Response", n = length(intersect(NoResp$ENTREZID,knowles5)))
Overlapk5 <- merge(Overlapk5,NRk5, all =TRUE)[,1:2]


length(intersect(intersect(knowles4, knowles5), NoResp$ENTREZID))

drug_palc <- c("#8B006D","#DF707E","#F1B72B", "#3386DD","#707031","#41B333")

ggplot(Overlapk4, aes(x=id, y=n))+
  geom_col(postition= "fill",aes(fill=id))+
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


#make as a set



knover45 <- list(Overlapk5, Overlapk4)
names(knover45) <- c("Knowles 4", "Knowles 5")
knover45 <- bind_rows(list("all eQTLs"=Overlapk4,"all response eQTLs"=Overlapk5), .id ="Overlap")
knover45$counts <- geneset_length <- rep(c(length(sigVDA24$ENTREZID),length(sigVDX24$ENTREZID),length(sigVEP24$ENTREZID), length(sigVMT24$ENTREZID),length(NoResp$ENTREZID)),2)


knover45 %>%
  mutate(percent=n/counts) %>%
  # ggplot(., aes(counts, y=n))+
  ggplot(., aes(x=Overlap, y=percent, fill = id))+
  geom_col(position = "fill")+
  geom_text(aes(label =  sprintf("%0.4f", round(percent, digits = 4))), position = position_stack(vjust=0.5))+
  scale_color_brewer(palette = "Dark2",guide = "none")+
  scale_fill_manual(values=drug_palc[c(1:4,6)])+
  theme_bw()+

  ylab("Percent")+
  xlab("")+
  ggtitle("Overlaps with Knowles data sets")+
  theme(plot.title = element_text(size=18,hjust = 0.5),
        axis.title = element_text(size = 15, color = "black"),
        axis.ticks = element_line(linewidth = 1.5),
        axis.line = element_line(linewidth = 1.5),
        axis.text.x = element_text(size = 12, color = "black", angle = 0),
        strip.text.x = element_text(size = 15, color = "black", face = "bold"))


