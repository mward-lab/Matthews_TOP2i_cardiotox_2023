# Cormotif
library(edgeR)
library(Cormotif)
library(RColorBrewer)
library(tidyverse)
library(BiocParallel)
## read in count file##
design <- read.csv("data/data_outline.txt", row.names = 1)
 # <- read.csv("data/cpmnorm_counts.csv",row.names = 1)
mymatrix <- readRDS("data/filtermatrix_x.RDS")#should be 14084
x_counts <- mymatrix$counts



indv <- as.factor(rep(c(1,2,3,4,5,6), c(12,12,12,12,12,12)))
time <- rep((rep(c("3h", "24h"), c(6,6))), 6)
time <- ordered(time, levels =c("3h", "24h"))
group <- as.factor(rep((c("1","2","3","4","5","6","7","8","9","10","11","12")),6))
drug <- rep(c("Daunorubicin","Doxorubicin","Epirubicin","Mitoxantrone","Trastuzumab", "Vehicle"),12)
group1 <- interaction(drug,time)
label <- (interaction(substring(drug, 0, 2), indv, time))
colnames(x_counts) <- label

# All data ----------------------------------------------------------------


group_fac <- group1
groupid <- as.numeric(group_fac)

compid <- data.frame(c1= c(1,2,3,4,5,7,8,9,10,11), c2 = c( 6,6,6,6,6,12,12,12,12,12))

y_TMM_cpm <- cpm(x_counts, log = TRUE)
#y_TMM_cpm <- filcpm_matrix
colnames(y_TMM_cpm) <- label
y_TMM_cpm
set.seed(12345)

##after execution, was saved to the RDS
# cormotif_initial <- cormotiffit(exprs = y_TMM_cpm,
                             # groupid = groupid,
                             # compid = compid,
                             # K=1:8, max.iter = 500,runtype="logCPM")
# cormotif_initialX <- cormotiffit(exprs = y_TMM_cpm,
#                                 groupid = groupid,
#                                 compid = compid,
#                                 K=5, max.iter = 400, runtype="logCPM")

# saveRDS(cormotif_initial,"data/cormotif_initialall.RDS")##saved so Ido not have to run every time
cormotif_initial <- readRDS("data/cormotif_initialall.RDS")##getting the final data back in

plotIC(cormotif_initial)

plotMotif(cormotif_initial)
# plot.new()
# legend("topright", legend = reverse(c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9", "1"), fill = my_barcolor)) #gray(seq(from = 1, to = 0, by = -0.1))[1:2])

gene_prob_tran <- cormotif_initial$bestmotif$p.post
rownames(gene_prob_tran) <- rownames(y_TMM_cpm)
dim(gene_prob_tran)

motif_prob <- cormotif_initial$bestmotif$clustlike
rownames(motif_prob) <- rownames(y_TMM_cpm)


clust1 <- motif_prob %>% ##NR
  as.data.frame() %>%
  filter(V1>0.5) %>%
  rownames
clust2 <- motif_prob %>% ##TI
  as.data.frame() %>%
  filter(V2>0.5) %>%
  rownames
clust3 <- motif_prob %>% ##ER
  as.data.frame() %>%
  filter(V3>0.5) %>%
  rownames
clust4 <- motif_prob %>% ##LR
  as.data.frame() %>%
  filter(V4>0.5) %>%
  rownames
intersect(clust3,clust4)



cormotif_initial$bestmotif$motif.prior
#[1] 0.52612387 0.04184878 0.03463675 0.39739060
#saveRDS(gene_prob_tran,"data/gene_probabilityk5.RDS")
##old name- non response
old_clust1  <- rownames(gene_prob_tran[(gene_prob_tran[,1] <0.45 &
                                                   gene_prob_tran[,2] <0.45 &
                                                   gene_prob_tran[,3] <0.45 &
                                                   gene_prob_tran[,4] <0.45 &
                                                   gene_prob_tran[,5] <0.45 &
                                                   gene_prob_tran[,6] <0.45 &
                                                   gene_prob_tran[,7] <0.45 &
                                                   gene_prob_tran[,8] <0.45 &
                                                   gene_prob_tran[,9] <0.45 &
                                                   gene_prob_tran[,10] <0.45),])
length(intersect(old_clust1,rownames(clust1)))




length((old_clust1)) ##7362
#sample(old_clust1,4)
#"55686"     "100128682" "11190"     "58472"
#"150786" "54856"  "145165" "23095"
#write_csv(as.data.frame(nonresponse_cluster), "data/cormotif_NRset.txt")

# TI_respset --------------------------------------------------------------
# all_response <-rownames(gene_prob_tran)
# response_set <- setdiff(all_response,nonresponse_cluster)
TI_respint1  <- rownames(gene_prob_tran[(gene_prob_tran[,1]> 0.10 &
                                         gene_prob_tran[,2] > 0.10 &
                                         gene_prob_tran[,3] > 0.10 &
                                         gene_prob_tran[,4] > 0.10 &
                                         # gene_prob_tran[,5] < 0.10 &
                                           gene_prob_tran[,6] > 0.10&
                                           gene_prob_tran[,7]> 0.10 &
                                           gene_prob_tran[,8]> 0.10 &
                                           gene_prob_tran[,9]> 0.10),])# &
                                           # gene_prob_tran[,10]< 0.10),])
length(intersect(TI_respint1,clust2))#251
 length(TI_respint1) #432
 length(all_response)

 length(intersect(TI_respint1,ER_respint1))

# storetem <-  gene_prob_tran[rownames(gene_prob_tran) %in% response_set,]
 #116
 #write_csv(as.data.frame(TI_respint1), "data/cormotif_TI_cluster.txt")
 sample(TI_respint1,4)
 ##[1] "90417"     "105373311" "2959"      "27145"
 #"51278" "51043" "26152" "79832"
# ER_resp set ------------------------------------------------------

# cormotif_initial$bestmotif$motif.prior
# cormotif_initial$bestmotif$motif.q
 ER_respint1  <- rownames(gene_prob_tran[(gene_prob_tran[,1]>0.55 &
                                           gene_prob_tran[,2] >0.55 &
                                           gene_prob_tran[,3] >0.55 &
                                           gene_prob_tran[,4] >0.25&
                                           gene_prob_tran[,5] <0.90 &
                                           gene_prob_tran[,6] <0.9 &
                                           gene_prob_tran[,7]<0.9&
                                           gene_prob_tran[,8]<0.9 &
                                          gene_prob_tran[,9]<0.9 &
                                           gene_prob_tran[,10]<0.90),])
                                           #
 length(unique(ER_respint1))
#481


#write_csv(as.data.frame(ER_respint1), "data/cormotif_ER_cluster.txt")

sample(ER_respint1,4)
#"54434"  "23506"  "51058"  "283219"
length(intersect(clust3,ER_respint1))#414
# LR_respset --------------------------------------------------------------

 LR_respint1  <- rownames(gene_prob_tran[(gene_prob_tran[,1] <0.970 &
                                            gene_prob_tran[,2] <0.97 &
                                            gene_prob_tran[,3] <0.97 &
                                            gene_prob_tran[,4] <0.97 &
                                            gene_prob_tran[,5] <0.9&
                                            gene_prob_tran[,6] >0.55 &
                                            gene_prob_tran[,7] >0.55 &
                                            gene_prob_tran[,8] >0.55 &
                                            gene_prob_tran[,9] >0.05 &
                                            gene_prob_tran[,10] <0.9),])

 length(LR_respint1)
#4850
 length(intersect(clust4,LR_respint1))#4675
# write_csv(as.data.frame(LR_respint1), "data/cormotif_LR_cluster.txt")
 sample(LR_respint1,4)
# "8915"   "1756"   "92335"  "146754"

# list summary ------------------------------------------------------------
filter(clust1 %in% )
length(nonresponse_cluster)
 #motif is 7362  set of gene with less than 50% chance of being DE

 length(TI_respint1)
 #432

 length(ER_respint1)
 #481
 #intersection of all early TOP2Bi genes that have greater than 50% chance of being DE

 length(LR_respint1)
 #4850
 #intersection of all TOP2Bi genes that have greater than 50% chance of being DE

 motif_TI <- clust2

 motif_ER <- clust3

 motif_LR <- clust4

 motif_NR <- clust1

 DEG_cormotif <- list(motif_NR,motif_TI,motif_LR,motif_ER)

 names(DEG_cormotif) <- c('motif_NR','motif_TI','motif_LR','motif_ER')

 saveRDS(DEG_cormotif, "data/DEG_cormotif.RDS")


