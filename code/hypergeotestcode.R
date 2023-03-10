#BiocManager::install("GeneOverlap")

library(fs)
library(readr)
library(tidyverse)
library(boot)
library(VennDiagram)
library(stats)
library(stats4)
#detach(package:dr4pl, unload = TRUE)

# BiocManager::install("montilab/hypeR", version="devel")
#
#
# library(hyperR)


# Data loading ------------------------------------------------------------
mymatrix <- readRDS('data/mymatrix.RDS')
allgenes <- mymatrix$genes[1]
# write_csv(allgenes, "data/allgenes.txt")

#load the files:
file.names <- list.files(path = "data/",
                         pattern = "sigV*", include.dir=FALSE,ignore.case = TRUE,recursive=FALSE,
                         full.names = TRUE)
file.names <- file.names[is_file(file.names)]
all_genes <- read_csv("data/allgenes.txt",show_col_types = FALSE)###complete
filenameonly <- read_csv("data/filenameonly.txt",show_col_types = FALSE)

for (k in 1:length(file.names)){

  assign(paste0(filenameonly$x[k]) , read.csv(file.names[k]))
}
backGL <- read_csv("data/backGL.txt",
                   col_types = cols(...1 = col_skip()))
backGL <- backGL %>% rename("ENTREZID"=1) %>% select(ENTREZID)
##rename the columns to the previous names
count(backGL) #14823
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
chrom_reg_Seoane <- read_csv(file = "data/Seonane2019supp1.txt",col_types = cols(...1 = col_skip()))

Seoane_2019 <- chrom_reg_Seoane[,2]
names(Seoane_2019) <- "ENTREZID"
total24 <-c(sigVDA24,sigVDX24,sigVEP24,sigVMT24)
AC_24list <- c(sigVDA24,sigVDX24,sigVEP24)
# intersections files ----------------------------------------------------------
list24totvenn <- get.venn.partitions(total24)
ACunion<- union(sigVDA24,sigVDX24)
ACunion <- union(ACunion,sigVEP24)
DDEresp <- list24totvenn$..values..[[9]]  #(4400)
DDEMresp <- list24totvenn$..values..[[1]]
Dxresp <- sigVDX24$ENTREZID
DXsprespon <- list24totvenn$..values..[[14]]

top2bi <-  read_csv("data/response_cluster24h.csv",
                    col_types = cols(...1 = col_skip()))

top2bi <- top2bi %>% rename("ENTREZID"=1)
length(top2bi$ENTREZID)



# functions ---------------------------------------------------------------

##will take a LONG TIME!
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


data_total24 <- overlap_significance(genes_all=backGL[[1]], gene_sets =total24, 5000)
# data_top2bi <- overlap_significance(genes_all=backGL[[1]], gene_sets =top2bi[[1]], 10000)
# data_DDEresp24<- overlap_significance(genes_all=backGL[[1]], gene_sets =DDEresp, 10000)
# data_DXresp24<- overlap_significance(genes_all=backGL[[1]], gene_sets =Dxresp, 10000)
# data_DX_only_resp24<- overlap_significance(genes_all=backGL[[1]], gene_sets =DXsprespon, 10000)
AC_24overlap <- overlap_significance(genes_all=backGL[[1]], gene_sets =AC_24list, 10000)

# plotting data_total24 ---------------------------------------------------

data.frame(test= data_total24$simulated_values) %>% #, Top2Bi_overlap=hyperDDEM24, Plot_Motif_set=hyperTop2bi, AC_24_overlap=hyperDDE24) %>%
  pivot_longer(everything(), names_to= "GeneList", values_to="Overlap") %>%
  ggplot(aes(x=Overlap))+
  geom_histogram(binwidth = 1)+
  annotate("text", x= 500, y=30, label =paste0("p =",data_total24$pval))+
  geom_vline(xintercept=data_total24$observed, lty=2, color="grey", size=1) +
  #geom_vline(xintercept=2201, lty=2, color="grey", size=1)+
  theme_bw() +
  theme(text=element_text(size=12), legend.position="none") +
  facet_grid(GeneList~.) +
  #scale_fill_manual(values=c("dodgerblue", "seagreen")) +
  ylab("Count")+
  ggtitle("Probabliity of overlap of all DEG at 24 hours, n=5000")


data.frame(test= AC_24overlap$simulated_values) %>% #, Top2Bi_overlap=hyperDDEM24, Plot_Motif_set=hyperTop2bi, AC_24_overlap=hyperDDE24) %>%
  pivot_longer(everything(), names_to= "GeneList", values_to="Overlap") %>%
  ggplot(aes(x=Overlap))+
  geom_histogram(binwidth = 1)+
  annotate("text", x= 3000, y=50, label =paste0("p =",AC_24overlap$pval))+
  geom_vline(xintercept=AC_24overlap$observed, lty=2, color="grey", size=1) +
  #geom_vline(xintercept=2201, lty=2, color="grey", size=1)+
  theme_bw() +
  theme(text=element_text(size=12), legend.position="none") +
  facet_grid(GeneList~.) +
  #scale_fill_manual(values=c("dodgerblue", "seagreen")) +
  ylab("Count")+
  ggtitle("Probabliity of overlap of all AC DEG at 24 hours n=10,000")





























#
# hyper <- rhyper(
#   nn=length(all_genes), m=length(sets[[1]]),
#   n=length(all_genes) - length(sets[[1]]),
#   k=length(sets[[2]])
# )
#
# simulated <- map_dbl(seq_len(10000), function(x) {
#   sim <- map(lengths(sets), ~sample(all_genes, .x))
#   sim <- length(reduce(sim, intersect))
#   return(sim)
# })
# data.frame(Simulated=simulated, Hypergeometric=hyper) %>%
#   pivot_longer(everything(), names_to="Method", values_to="Overlap") %>%
#   ggplot(aes(x=Overlap, fill=Method)) +
#   geom_histogram(binwidth=1) +
#   geom_vline(xintercept=observed, lty=2, color="grey", size=1) +
#   theme_bw() +
#   theme(text=element_text(size=12), legend.position="none") +
#   facet_grid(Method~.) +
#   scale_fill_manual(values=c("dodgerblue", "seagreen")) +
#   ylab("Count")
# sim_pval <- (sum(simulated >= observed) + 1) / (10000 + 1)
# hyper_pval <- phyper(
#   q=observed, m=setsofgenes,
#   n=length(all_genes) - length(sets[[1]]),
#   k=length(sets[[2]]), lower.tail=FALSE
# )
# k=
# # Define your gene lists
# all_genes <- allgenes[[1]]
# subset_genes <- top2bi[[1]]
#
# # Define your background
# total_genes <- length(all_genes)
#
# # Set the parameters for the hypergeometric test
# k <- length(subset_genes)   # number of successes in the sample
# M <- total_genes            # population size
# n <- k                      # sample size
#
# # Calculate the p-value using the phyper() function
# pval <- 1 - phyper(k-1, M-k, total_genes-M, n, lower.tail=FALSE)
#
# # Print the p-value
# print(pval)
#
# # Define the parameters
# x <- 4400      # Number of successes in the sample
# m <- 8761  # Total number of objects in the population
# n <- 6793      # Sample size (total -black balls)
# k <- 2201:4400    # Number of objects classified as successes
# bb <- 6793  #number of black balls
# # Calculate the probability mass function
# pmf <- sum(dhyper(x, m, k, n))
#
# # Print the result
# print(pmf)
# #The dhyper() function requires the user to specify the number of successes in
# #the sample (x), the total number of objects in the population (m), the sample
# #size (n), and the number of objects in the population that are classified as
# #successes (k). The function then calculates the probability of observing
# #exactly x successes in a sample of size n drawn from a population of size m
# #with k successes in the population.
#
# random <- rhyper(nn=10000,m,n,k)
#
# median(random)
#
#
#
# #observed <- length(DDEresp)
# data.frame(test= random) %>% #, Top2Bi_overlap=hyperDDEM24, Plot_Motif_set=hyperTop2bi, AC_24_overlap=hyperDDE24) %>%
#   pivot_longer(everything(), names_to= "GeneList", values_to="Overlap") %>%
#                  ggplot(aes(x=Overlap))+
#                     geom_histogram(binwidth = 1)+
#   annotate("text", x=3000, y=20, label =paste0("p =",pmf))+
#                     geom_vline(xintercept=4400, lty=2, color="grey", size=1) +
#   geom_vline(xintercept=2201, lty=2, color="grey", size=1)+
#                     theme_bw() +
#                     theme(text=element_text(size=12), legend.position="none") +
#                     facet_grid(GeneList~.) +
#                     #scale_fill_manual(values=c("dodgerblue", "seagreen")) +
#                     ylab("Count")+
#   ggtitle(paste0("Random distribution choosing ",k," balls from ",m," balls\n containing ",bb," black balls "))
# qhyper(0.5:1, m,n,k)
# phyper(q=4400, m,n,k)
#
