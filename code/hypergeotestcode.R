#BiocManager::install("GeneOverlap")


library(readr)
library(tidyverse)
library(boot)
library(VennDiagram)
#detach(package:dr4pl, unload = TRUE)

BiocManager::install("montilab/hypeR", version="devel")


library(hyperR)


#load the files:
file.names <- list.files(path = "data/",
                         pattern = "sig*", ignore.case = TRUE,
                         full.names = TRUE)

filenameonly <- read_csv("data/filenameonly.txt")


for (k in 1:length(file.names)){

  assign(paste0(filenameonly$x[k]) , read.csv(file.names[k]))
}
backGL <- read_csv("data/backGL.txt",
                   col_types = cols(...1 = col_skip()))
backGL <- backGL %>% rename("ENTREZID"=1) %>% select(ENTREZID)
##rename the columns to the previous names
count(backGL)
sigVDA24<- sigVDA24 %>% rename("ENTREZID"=1) %>% select(ENTREZID)
sigVDX24<- sigVDX24 %>% rename("ENTREZID"=1) %>% select(ENTREZID)
sigVEP24<- sigVEP24%>% rename("ENTREZID"=1) %>% select(ENTREZID)
sigVMT24<- sigVMT24%>% rename("ENTREZID"=1) %>% select(ENTREZID)
sigVTR24<- sigVTR24%>% rename("ENTREZID"=1) %>% select(ENTREZID)
sigVDA3<- sigVDA3%>% rename("ENTREZID"=1) %>% select(ENTREZID)
sigVDX3<- sigVDX3%>% rename("ENTREZID"=1) %>% select(ENTREZID)
sigVEP3<- sigVEP3%>% rename("ENTREZID"=1) %>% select(ENTREZID)
sigVMT3<- sigVMT3%>% rename("ENTREZID"=1) %>% select(ENTREZID)
sigVTR3 <-sigVTR3%>% rename("ENTREZID"=1) %>% select(ENTREZID)


total24 <-c(sigVDA24,sigVDX24,sigVEP24,sigVMT24)
AC_24list <- c(sigVDA24,sigVDX24,sigVEP24)
# intersections files ----------------------------------------------------------
list24totvenn <- get.venn.partitions(total24)


DDEresp <- list24totvenn$..values..[[9]]  #(4400)
DDEMresp <- list24totvenn$..values..[[1]]
Dxresp <- sigVDX24$ENTREZID
DXsprespon <- list24totvenn$..values..[[14]]

top2bi <-  read_csv("data/response_cluster24h.csv",
                    col_types = cols(...1 = col_skip()))

top2bi <- top2bi %>% rename("ENTREZID"=1)
length(top2bi$ENTREZID)



# functions ---------------------------------------------------------------


overlap_significance <- function(genes_all, gene_sets, iterations) {
  observed <- length(reduce(gene_sets, intersect))
  simulated <- map_dbl(seq_len(iterations), function(x) {
    sim <- map(lengths(gene_sets), ~sample(genes_all, .x))
    sim <- length(reduce(sim, intersect))
    return(sim)
  })
  pval <- (sum(simulated >= observed) + 1) / (iterations + 1)
  return(list(pval=pval, simulated_values=simulated, observed=observed))
}

data_total24 <- overlap_significance(genes_all=backGL[[1]], gene_sets =total24, 10000)
data_top2bi <- overlap_significance(genes_all=backGL[[1]], gene_sets =top2bi[[1]], 10000)
data_DDEresp24<- overlap_significance(genes_all=backGL[[1]], gene_sets =DDEresp, 10000)
data_DXresp24<- overlap_significance(genes_all=backGL[[1]], gene_sets =Dxresp, 10000)
data_DX_only_resp24<- overlap_significance(genes_all=backGL[[1]], gene_sets =DXsprespon, 10000)


observed <- length(reduce(ac_overlap24, intersect))

observed

ggplot(data_total24, aes(x=simulated_vlues))+geom_histogram()
all_genes <- backGL#sprintf("ENSG%08d", seq_len(10000))
common_genes <-DDEresp

head(all_genes, 5)

sets <- map(seq_len(2), function(x) {
  c(common_genes, sample(all_genes, sample(seq(200, 500), 1)))
})

 lengths(sets)

observed <- length(reduce(sets, intersect))

hyper <- rhyper(
  nn=length(all_genes), m=length(sets[[1]]),
  n=length(all_genes) - length(sets[[1]]),
  k=length(sets[[2]])
)

simulated <- map_dbl(seq_len(10000), function(x) {
  sim <- map(lengths(sets), ~sample(all_genes, .x))
  sim <- length(reduce(sim, intersect))
  return(sim)
})
data.frame(Simulated=simulated, Hypergeometric=hyper) %>%
  pivot_longer(everything(), names_to="Method", values_to="Overlap") %>%
  ggplot(aes(x=Overlap, fill=Method)) +
  geom_histogram(binwidth=1) +
  geom_vline(xintercept=observed, lty=2, color="grey", size=1) +
  theme_bw() +
  theme(text=element_text(size=12), legend.position="none") +
  facet_grid(Method~.) +
  scale_fill_manual(values=c("dodgerblue", "seagreen")) +
  ylab("Count")
sim_pval <- (sum(simulated >= observed) + 1) / (10000 + 1)
hyper_pval <- phyper(
  q=observed, m=length(sets[[1]]),
  n=length(all_genes) - length(sets[[1]]),
  k=length(sets[[2]]), lower.tail=FALSE
)
sets <-

hyperDDE24 <- rhyper(
  nn=length(backGL[[1]]), m=length(DDEresp),
  n=length(backGL[[1]]) - length(DDEresp),
  k=10000)
hyperDDE24_pval <- phyper(
  q=(length(DDEresp)-1), m=length(DDEresp),
  n=length(backGL[[1]]) - length(DDEresp),
  k=length(DDEresp), lower.tail=FALSE)
hyperDDE24_pval
quantile(hyperDDE24)
contable <- matrix(c(
  dg=length(signature)-length(DDEresp),
  dr=length(DDEresp),
  ng=backpop-length(signature)-length(pathway)+length(DDEresp),
  nr=length(pathway)-length(DDEresp)),2,2,dimnames=list(c("GREEN","RED"),c("DRAWN","not DRAWN")))

in_common24 <-c(sigVDA24$ENTREZID,sigVDX24$ENTREZID,sigVEP24$ENTREZID,sigVMT24$ENTREZID)
incommon_24 <- (unique(in_common24))



hyperTop2bi <- rhyper(
  nn=length(backGL[[1]]), m=length(top2bi[[1]]),
  n=length(backGL[[1]]) - length(top2bi[[1]]),
  k=10000)

# hyperDx24all <- rhyper(
#   nn=length(incommon_24), m=length(Dxresp),
#   n=length(incommon_24) - length(Dxresp),
#   k=10000)

hyperDDEM24 <- rhyper(
  nn=length(DDEMresp), m=length(DDEMresp),
  n=length(incommon_24) - length(DDEMresp),
  k=10000)

length(Dxresp)

observed <- length(DDEresp)
data.frame(test= hyperDDE24) %>% #, Top2Bi_overlap=hyperDDEM24, Plot_Motif_set=hyperTop2bi, AC_24_overlap=hyperDDE24) %>%
  pivot_longer(everything(), names_to= "GeneList", values_to="Overlap") %>%
                 ggplot(aes(x=Overlap))+
                    geom_histogram(binwidth = 1)+
                    #geom_vline(xintercept=observed, lty=2, color="grey", size=1) +
                    theme_bw() +
                    theme(text=element_text(size=12), legend.position="none") +
                    facet_grid(GeneList~.) +
                    #scale_fill_manual(values=c("dodgerblue", "seagreen")) +
                    ylab("Count")

