# Cormotif
library(edgeR)
library(Cormotif)
library(RColorBrewer)
## read in count file##
design <- read.csv("data/data_outline.txt", row.names = 1)
 # <- read.csv("data/cpmnorm_counts.csv",row.names = 1)
mymatrix <- readRDS("data/filtermatrix_x.RDS")#should be 14084
x_counts <- mymatrix$counts
colSums(x_counts)


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
y_TMM_cpm <- filcpm_matrix
colnames(y_TMM_cpm) <- label
y_TMM_cpm
set.seed(12345)
cormotif_initial <- cormotiffit(exprs = y_TMM_cpm,
                             groupid = groupid,
                             compid = compid_tran,
                             K=1:8, max.iter = 500)
cormotif_5 <- cormotiffit(exprs = y_TMM_cpm,
                                groupid = groupid,
                                compid = compid_tran,
                                K=5, max.iter = 400)

saveRDS(cormotif_initial,"data/cormotif_initialall.RDS")##saved so Ido not have to run everytime
#cormotif_initial <- readRDS("data/cormotif_initialK5.RDS")

plotIC(cormotif_initial)



plotMotif(cormotif_initial)
# plot.new()
# legend("topright", legend = reverse(c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9", "1"), fill = my_barcolor)) #gray(seq(from = 1, to = 0, by = -0.1))[1:2])

gene_prob_tran <- cormotif_initial$bestmotif$p.post
rownames(gene_prob_tran) <- rownames(y_TMM_cpm)
dim(gene_prob_tran)

cormotif_initial$bestmotif$motif.prior
#[1] 0.52607840 0.04183696 0.03462484 0.00010635 0.39735345
#saveRDS(gene_prob_tran,"data/gene_probabilityk5.RDS")

nonresponse_cluster  <- rownames(gene_prob_tran[(gene_prob_tran[,1] <0.45 &
                                                   gene_prob_tran[,2] <0.45 &
                                                   gene_prob_tran[,3] <0.45 &
                                                   gene_prob_tran[,4] <0.45 &
                                                   gene_prob_tran[,5] <0.45 &
                                                   gene_prob_tran[,6] <0.45 &
                                                   gene_prob_tran[,7] <0.45 &
                                                   gene_prob_tran[,8] <0.45 &
                                                   gene_prob_tran[,9] <0.45 &
                                                   gene_prob_tran[,10] <0.45),])

length((nonresponse_cluster)) ##7362
sample(nonresponse_cluster,4)
#"55686"     "100128682" "11190"     "58472"
#"150786" "54856"  "145165" "23095"
write_csv(as.data.frame(nonresponse_cluster), "data/cormotif_NRset.txt")

# TI_respset --------------------------------------------------------------
all_response <-rownames(gene_prob_tran)
response_set <- setdiff(all_response,nonresponse_cluster)
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

 length(TI_respint1)
 length(all_response)

 length(intersect(TI_respint1,response_set))

storetem <-  gene_prob_tran[rownames(gene_prob_tran) %in% response_set,]
 #116
 write_csv(as.data.frame(TI_respint1), "data/cormotif_TI_cluster.txt")
 sample(TI_respint1,4)
 ##[1] "90417"     "105373311" "2959"      "27145"
 #"51278" "51043" "26152" "79832"
# ER_resp set ------------------------------------------------------

cormotif_initial$bestmotif$motif.prior
cormotif_initial$bestmotif$motif.q
 ER_respint1  <- rownames(gene_prob_tran[(gene_prob_tran[,1]>0.55 &
                                           gene_prob_tran[,2] >0.55 &
                                           gene_prob_tran[,3] >0.55 &
                                           gene_prob_tran[,4] >0.35&
                                           gene_prob_tran[,5] <0.90 &
                                           gene_prob_tran[,6] <0.9 &
                                           gene_prob_tran[,7]<0.9&
                                           gene_prob_tran[,8]<0.9 &
                                          gene_prob_tran[,9]<0.9 &
                                           gene_prob_tran[,10]<0.90),])
                                           #
 length(unique(ER_respint1))
#422


write_csv(as.data.frame(ER_respint1), "data/cormotif_ER_cluster.txt")

sample(ER_respint1,4)
#"54434"  "23506"  "51058"  "283219"
length(intersect(TI_respint1,ER_respint1))
# LR_respset --------------------------------------------------------------

 LR_respint1  <- rownames(gene_prob_tran[(gene_prob_tran[,1] <0.50 &
                                            gene_prob_tran[,2] <0.50 &
                                            gene_prob_tran[,3] <0.50 &
                                            gene_prob_tran[,4] <0.50 &
                                            gene_prob_tran[,5] <0.5&
                                            gene_prob_tran[,6] >0.10 &
                                            gene_prob_tran[,7] >0.10 &
                                            gene_prob_tran[,8] >0.10 &
                                            gene_prob_tran[,9] >0.10 &
                                            gene_prob_tran[,10] <0.5),])

 length(unique(LR_respint1))
 #2929
 length(intersect(ER_respint1,LR_respint1))
 write_csv(as.data.frame(LR_respint1), "data/cormotif_LR_cluster.txt")
 sample(LR_respint1,4)
# "8915"   "1756"   "92335"  "146754"
# motif2 ------------------------------------------------------------------

 motif2_unknown  <- rownames(gene_prob_tran[(gene_prob_tran[,1] >0.06 &
                                            gene_prob_tran[,2] >0.06 &
                                            gene_prob_tran[,3] >0.06 &
                                            gene_prob_tran[,4] >0.04 &
                                            gene_prob_tran[,5] >0.02 &
                                            gene_prob_tran[,6] >0.06 &
                                            gene_prob_tran[,7] >0.06 &
                                            gene_prob_tran[,8] >0.06 &
                                            gene_prob_tran[,9] >0.04 &
                                            gene_prob_tran[,10] >0.02),])

 length(unique(LR_respint1))
 write_csv(as.data.frame(LR_respint1), "data/cormotif_LR_respint.txt")
#22822,4356


# list summary ------------------------------------------------------------

nonresponse_cluster
 #motif is 8705  set of gene with less than 50% chance of being DE

 TI_respint1
 #116
 #intersection of all TOP2Bi genes that have greater than 50% chance of being DE and less than in Tras
 ER_respint1
 #424
 #intersection of all early TOP2Bi genes that have greater than 50% chance of being DE

 LR_respint1
 #1409
 #intersection of all TOP2Bi genes that have greater than 50% chance of being DE

 motif_TI <- TI_respint1
#116
 motif_ER <- ER_respint1
#452
 motif_LR <- LR_respint1
#1515
 motif_NR <- nonresponse_cluster
#8846
 DEG_cormotif <- list(motif_NR,motif_TI,motif_LR,motif_ER)

 names(DEG_cormotif) <- c('motif_NR','motif_TI','motif_LR','motif_ER')

 saveRDS(DEG_cormotif, "data/DEG_cormotif.RDS")




# limmafit counts ---------------------------------------------------------

testing(y_TMM_cpm,groupid,compid)

# Load gene counts data
data(countdata)

# Cluster the data
clust_result <- cormotif(countdata)

# Extract cluster assignments
cluster_assignments <- clusterAssign(clust_result)

# Extract gene names for each cluster
gene_names <- lapply(cluster_assignments, function(x) rownames(countdata)[x])







# old ideas-not in use as of march 20 --------------------------------------




#generesponse tran contains all probabilities

 #use cormotif_initial$bestmotif$motif.q for  probabilities to subset the columns by

q5geneset <- rownames(gene_prob_tran[(gene_prob_tran[,1]>q5[1]<0.5&
                            gene_prob_tran[,2] >q5[2]&
                            gene_prob_tran[,3] >q5[3]&
                            gene_prob_tran[,4] >q5[4]&
                            gene_prob_tran[,5] >q5[5]&
                            gene_prob_tran[,6] >q5[6]&
                            gene_prob_tran[,7] >q5[7]&
                            gene_prob_tran[,8] >q5[8]&
                            gene_prob_tran[,9] >q5[9]&
                            gene_prob_tran[,10]>q5[10]),])
 q4geneset <- rownames(gene_prob_tran[(gene_prob_tran[,1]> q2[1]&
                                         gene_prob_tran[,2] >q2[2]&
                                         gene_prob_tran[,3] >q2[3]&
                                         gene_prob_tran[,4] >q2[4]&
                                         gene_prob_tran[,5] >q2[5]&
                                         gene_prob_tran[,6] >q2[6]&
                                         gene_prob_tran[,7] >q2[7]&
                                         gene_prob_tran[,8] >q2[8]&
                                         gene_prob_tran[,9] >q2[9]&
                                         gene_prob_tran[,10]>
                                         q3[10]),])
 q3geneset <- rownames(gene_prob_tran[(gene_prob_tran[,1]> q2[1]&
                                         gene_prob_tran[,2] >q2[2]&
                                         gene_prob_tran[,3] >q2[3]&
                                         gene_prob_tran[,4] >q2[4]&
                                         gene_prob_tran[,5] >q2[5]&
                                         gene_prob_tran[,6] >q2[6]&
                                         gene_prob_tran[,7] >q2[7]&
                                         gene_prob_tran[,8] >q2[8]&
                                         gene_prob_tran[,9] >q2[9]&
                                         gene_prob_tran[,10]>
                                         q3[10]),])

 q2geneset <- rownames(gene_prob_tran[(gene_prob_tran[,1]> q2[1]&
                                         gene_prob_tran[,2] >q2[2]&
                                         gene_prob_tran[,3] >q2[3]&
                                         gene_prob_tran[,4] >q2[4]&
                                         gene_prob_tran[,5] <q2[5]&
                                         gene_prob_tran[,6] >q2[6]&
                                         gene_prob_tran[,7] >q2[7]&
                                         gene_prob_tran[,8] >q2[8]&
                                         gene_prob_tran[,9] >q2[9]&
                                         gene_prob_tran[,10]<
                                         q3[10]),])
 q1geneset <- rownames(gene_prob_tran[(gene_prob_tran[,1]< q1[1]&
                                         gene_prob_tran[,2] <q1[2]&
                                         gene_prob_tran[,3] <q1[3]&
                                         gene_prob_tran[,4] <q1[4]&
                                        gene_prob_tran[,5] <q2[5]&
                                         gene_prob_tran[,6] <q1[6]&
                                         gene_prob_tran[,7] <q1[7]&
                                         gene_prob_tran[,8] <q1[8]&
                                         gene_prob_tran[,9] <q1[9]&
                                         gene_prob_tran[,10]<
                                         q2[10]),])
 q3hrneset <- rownames(gene_prob_tran[(gene_prob_tran[,1]> 0.55[1]&
                                         gene_prob_tran[,2] >0.55[2]&
                                         gene_prob_tran[,3] >0.55[3]&
                                         gene_prob_tran[,4] >0.55[4]&
                                         gene_prob_tran[,5] <0.55[5]&
                                         gene_prob_tran[,6] >0.55[6]&
                                         gene_prob_tran[,7] >0.55[7]&
                                         gene_prob_tran[,8] >0.55[8]&
                                         gene_prob_tran[,9] >0.55[9]&
                                         gene_prob_tran[,10]<
                                         0.55[10]),])
testset <- (gene_prob_tran[(gene_prob_tran[,5]*(q1[5])<0.5),])
length(testset)
setdiff(q3geneset,q4geneset,)
 list(list_ofrank[which(list_ofrank == 1,arr.ind = TRUE)[,1],])
 row.names(df1[which(df1==17,arr.ind=T)[,1],] )
q1 <- c(0.0003493127, 0.0002418384, 0.0002650427, 0.0002566342, 0.0013186903 ,0.0003948442,
        0.0002322792 ,0.0002054453 ,0.0128938596 ,0.0013603579)
 q2 <-  c( 0.99741838, 0.97630177, 0.99704198, 0.98945155, 0.02508136, 0.99780119, 0.99730884, 0.99550571,
            0.27375237, 0.02282597)
 q3 <- c(0.997354673, 0.928833988, 0.996832588, 0.979228124, 0.028695405, 0.007432923, 0.003184927,
         0.003448530, 0.040165153, 0.022335571)
 q4 <- c(0.5296498, 0.4803906, 0.5079664, 0.4910399, 0.4948393, 0.5209333, 0.5007303, 0.4796137, 0.4951372,
         0.4953019)
 q5 <- c(0.0005805927, 0.0003467514, 0.0003414524, 0.0003861847, 0.0018299596, 0.9997636975,
          0.9997309701, 0.9993176883, 0.3496581051, 0.0020936440)
 cormotif_initial$bestmotif$motif.q[1,]
 cormotif_initial$bestmotif$motif.prior

 list_ofrank <- generank(cormotif_initial$bestmotif$p.post)

 ### none work so far do eyeballing
 # set1- no response_cluster24h 11190
 # set3- all top2bi 1275442
 # set4-late responsetop2bi only 126820
 # set5-early response top2bi only 27245
 #
 # geneexpressionsets <- cbind(sets=c('set1', 'set3','set4','set5'), ENTREZID = c(11190,1275442,126820,27245))
cormotif_initial$bestmotif$motif.prior


# Three hour patterns ------------------------------------------------------

#note: learned that I can just include the compid frames only to subset the right counts
group_fac <- group1
groupid <- as.numeric(group_fac)
h3groupid <- factor(rep(1:6,6),levels = c(1,2,3,4,5,6))
h3compid <- data.frame(c1= c(1,2,3,4,5), c2 = c( 6,6,6,6,6))

y_TMM_cpm <- cpm(x_counts, log = TRUE)


colnames(y_TMM_cpm) <- label
threehour <- y_TMM_cpm[,c(1:6,13:18,25:30,37:42,49:54,61:66)]
set.seed(12345)
cormotif_3h <- cormotiffit(exprs = threehour,
                               groupid = h3groupid,
                               compid = h3compid,
                               K=1:8, max.iter = 500)

plotIC(cormotif_3h)
#x_axis_labels(labels = c("3_Daun","3_Dox","3_Epi","3_Mito","3_Tras","24_Daun","24_Dox","24_Epi","24_Mito","24_Tras"), every_nth = 1, adj=1, srt =90, cex =0.4)
plotMotif(cormotif_3h)
saveRDS(cormotif_3h,"data/cormotif_3hk1-8.RDS")
gene_prob_tran3h <- cormotif_3h$bestmotif$p.post
rownames(gene_prob_tran3h) <- rownames(y_TMM_cpm)
dim(gene_prob_tran3h)
saveRDS(gene_prob_tran3h,"data/gene_prob_tran3h.RDS")
gene_prob_tran3h <- readRDS("data/gene_prob_tran3h.RDS")





NR_cluster3h  <- rownames(gene_prob_tran3h[(gene_prob_tran3h[,1] <0.5 &
                                              gene_prob_tran3h[,2] <0.5 &
                                              gene_prob_tran3h[,3] <0.5 &
                                              gene_prob_tran3h[,4] <0.5 &
                                              gene_prob_tran3h[,5] <0.5),])

length(NR_cluster3h)

#24 hour patterns -------------------------------------

compid_24h <- data.frame(c1= c(7,8,9,10,11), c2 = c( 12,12,12,12,12))

y_TMM_cpm <- cpm(x_counts, log = TRUE)


colnames(y_TMM_cpm) <- label
#threehour <- y_TMM_cpm[,c(1:6,13:18,25:30,37:42,49:54,61:66)]
set.seed(12345)
cormotif_24h2mot <- cormotiffit(exprs = y_TMM_cpm,
                           groupid = groupid,
                           compid = compid_24h,
                           K=1:5, max.iter = 500)
saveRDS(cormotif_24h2mot,file = "data/Cormotif_24_k1-5_raw.RDS")
cormotif_24h2mot <- readRDS("data/Cormotif_24_k1-5_raw.RDS")

#x_axis_labels(labels = c("3_Daun","3_Dox","3_Epi","3_Mito","3_Tras","24_Daun","24_Dox","24_Epi","24_Mito","24_Tras"), every_nth = 1, adj=1, srt =90, cex =0.4)
# plotMotif(cormotif_24h2mot)
#  gene_prob_tran24h <- cormotif_24h2mot$bestmotif$p.post
# rownames(gene_prob_tran24h) <- rownames(y_TMM_cpm)
# dim(gene_prob_tran24h)
# nonresponse_cluster24h <- rownames(gene_prob_tran24h[(gene_prob_tran24h[,1] <0.5 &
#                                                         gene_prob_tran24h[,2] <0.5 &
#                                                         gene_prob_tran24h[,3] <0.5 &
#                                                         gene_prob_tran24h[,4] <0.5&
#                                                         gene_prob_tran24h[,5] <0.5),])
# length(nonresponse_cluster24h)
# response_cluster24h <- rownames(gene_prob_tran24h[(gene_prob_tran24h[,1] > 0.5 &
#                                                      gene_prob_tran24h[,2] > 0.5 &
#                                                      gene_prob_tran24h[,3] > 0.5 &
#                                                      gene_prob_tran24h[,4] > 0.5
#                                                     ),])
#
# ACresponse_cluster24h <- rownames(gene_prob_tran24h[(gene_prob_tran24h[,1] >0.5 & gene_prob_tran24h[,2] >0.5 & gene_prob_tran24h[,3] >0.5),])
# # length(response_cluster24h)
# # length(ACresponse_cluster24h)
# # write.csv(response_cluster24h,"data/Top2biresp_cluster24h.csv")
# # response_cluster24h <- read_csv("data/Top2biresp_cluster24h.csv")
# # length(intersect(allresponse_cluster,response_cluster24h[[2]]))
# #
# #
# # write.csv(ACresponse_cluster24h,"data/resp_cluster24h.csv")
# # write.csv(nonresponse_cluster24h,"data/nonresponse_cluster24h.csv")
# # "Da.1.3h",  "Do.1.3h",  "Ep.1.3h",  "Mi.1.3h",  "Tr.1.3h",  "Ve.1.3h",  "Da.1.24h", "Do.1.24h", "Ep.1.24h", "Mi.1.24h", "Tr.1.24h", "Ve.1.24h",
# # "Da.2.3h",  "Do.2.3h" , "Ep.2.3h", "Mi.2.3h" , "Tr.2.3h",,  "Ve.2.3h",  "Da.2.24h", "Do.2.24h", "Ep.2.24h", "Mi.2.24h", "Tr.2.24h", "Ve.2.24h",
# # "Da.3.3h",  "Do.3.3h",  "Ep.3.3h" , "Mi.3.3h" , "Tr.3.3h",  "Ve.3.3h",  "Da.3.24h", "Do.3.24h", "Ep.3.24h", "Mi.3.24h", "Tr.3.24h", "Ve.3.24h",
# # "Da.4.3h", "Do.4.3h",  "Ep.4.3h" , "Mi.4.3h" , "Tr.4.3h" , "Ve.4.3h",  "Da.4.24h", "Do.4.24h", "Ep.4.24h",  "Mi.4.24h", "Tr.4.24h", "Ve.4.24h",
# # "Da.5.3h",  "Do.5.3h",  "Ep.5.3h",  "Mi.5.3h",  "Tr.5.3h" , "Ve.5.3h",  "Da.5.24h", "Do.5.24h", "Ep.5.24h", "Mi.5.24h", "Tr.5.24h", "Ve.5.24h",
# # "Da.6.3h",  "Do.6.3h",  "Ep.6.3h",  "Mi.6.3h",  "Tr.6.3h",  "Ve.6.3h",  "Da.6.24h", "Do.6.24h", "Ep.6.24h", "Mi.6.24h", "Tr.6.24h", "Ve.6.24h")) {
#
# ## compute posterior cluster membership
# check<- limmafit(y_TMM_cpm, groupid,compid)
# ## Fit limma model
# limmafit<-function(exprs,groupid,compid) {
#   compnum<-nrow(compid)
#   genenum<-nrow(exprs)
#   limmat<-matrix(0,genenum,compnum)
#   limmas2<-rep(0,compnum)
#   limmadf<-rep(0,compnum)
#   limmav0<-rep(0,compnum)
#   limmag1num<-rep(0,compnum)
#   limmag2num<-rep(0,compnum)
#
#   for(i in 1:compnum) {
#     selid1<-which(groupid == compid[i,1])
#     selid2<-which(groupid == compid[i,2])
#     eset<-new("ExpressionSet", exprs=cbind(exprs[,selid1],exprs[,selid2]))
#     g1num<-length(selid1)
#     g2num<-length(selid2)
#     designmat<-cbind(base=rep(1,(g1num+g2num)), delta=c(rep(0,g1num),rep(1,g2num)))
#     fit<-lmFit(eset,designmat)
#     fit<-eBayes(fit)
#     limmat[,i]<-fit$t[,2]
#     limmas2[i]<-fit$s2.prior
#     limmadf[i]<-fit$df.prior
#     limmav0[i]<-fit$var.prior[2]
#     limmag1num[i]<-g1num
#     limmag2num[i]<-g2num
#
#     # log odds
#     # w<-sqrt(1+fit$var.prior[2]/(1/g1num+1/g2num))
#     # log(0.99)+dt(fit$t[1,2],g1num+g2num-2+fit$df.prior,log=TRUE)-log(0.01)-dt(fit$t[1,2]/w, g1num+g2num-2+fit$df.prior, log=TRUE)+log(w)
#   }
#   limmacompnum<-nrow(compid)
#   result<-list(t=limmat, v0=limmav0, df0=limmadf, s20=limmas2, g1num=limmag1num, g2num=limmag2num,compnum=limmacompnum)
#   #result<-list(t=limmat, v0=limmav0, df0=limmadf, s20=limmas2, g1num=limmag1num, g2num=limmag2num)
# }
# limmafit_custom(y_TMM_cpm,annotation=anno)
# limmafit_custom<-function(counts, annotation, time = c("1", "2"),
#                           drug = c("Daunorubicin", "Doxorubicin" , "Epirubicin",   "Mitoxantrone", "Trastuzumab","Vehcile")){
#   limmat <- vector()
#   limmas2 <- vector()
#   limmadf <- vector()
#   limmav0 <- vector()
#   limmag1num <- vector()
#   limmag2num <- vector()
#   limmacompnum <- length(drug) * length(time)
#
#   for (tp in time) {
#     for (b in drug) {
#       # if (tp == "48" & b == "Staph") {
#       #   limmacompnum <- limmacompnum - 1
#       #   next
#       }
#       message("Testing treatment", b, " at timepoint ", tp)
#       counts_sub <- counts[, annotation$drug %in% c("none", b) &
#                              annotation$time == tp]
#       anno_sub <- annotation[annotation$drug %in% c("none", b) &
#                                annotation$time == tp, ]
#       anno_sub <- droplevels(anno_sub)
#       library("edgeR")
#       params <- paste(anno_sub$drug, anno_sub$time, sep = ".")
#       params <- factor(params)
#       params <- relevel(params, ref = grep("none", levels(params), value = TRUE))
#       design <- model.matrix(~0 + params )
#       colnames(design) <- levels(params)
#       y <- DGEList(counts_sub)
#       y <- calcNormFactors(y)
#       v <- voom(y, design)
#       corfit <- duplicateCorrelation(v, design, block = anno_sub$indv)
#       v <- voom(y, design, block = anno_sub$ind, correlation = corfit$consensus)
#       fit <- lmFit(v, design, block = anno_sub$ind, correlation = corfit$consensus)
#       fit2 <- contrasts.fit(fit, contrasts = c(-1, 1, 0))
#       fit2 <- eBayes(fit2)
#       limmat <- cbind(limmat, fit2$t)
#       limmas2 <-c(limmas2, fit2$s2.prior)
#       limmadf <- c(limmadf, fit2$df.prior)
#       limmav0 <- c(limmav0, fit2$var.prior)
#       limmag1num <- c(limmag1num, sum(fit2$design[, 1]))
#       limmag2num <- c(limmag2num, sum(fit2$design[, 2]))
#     }
#   }
#
#   return(list(t=limmat, v0=limmav0, df0=limmadf, s20=limmas2, g1num=limmag1num, g2num=limmag2num,compnum=limmacompnum))
# }
#
# ## Log-likelihood for moderated t under H0
# modt.f0.loglike<-function(x,df) {
#   a<-dt(x, df, log=TRUE)
#   result<-as.vector(a)
#   flag<-which(is.na(result)==TRUE)
#   result[flag]<-0
#   result
# }
#
# ## Log-likelihood for moderated t under H1
# ## param=c(df,g1num,g2num,v0)
# modt.f1.loglike<-function(x,param) {
#   df<-param[1]
#   g1num<-param[2]
#   g2num<-param[3]
#   v0<-param[4]
#   w<-sqrt(1+v0/(1/g1num+1/g2num))
#   dt(x/w, df, log=TRUE)-log(w)
#   a<-dt(x/w, df, log=TRUE)-log(w)
#   result<-as.vector(a)
#   flag<-which(is.na(result)==TRUE)
#   result[flag]<-0
#   result
# }
#
# ## Correlation Motif Fit
# cmfit<-function(x, type, K=5, tol=1e-3, max.iter=100) {
#   ## initialize
#   xrow<-nrow(x)
#   xcol<-ncol(x)
#   loglike0<-list()
#   loglike1<-list()
#   p<-rep(1,K)/K
#   q<-matrix(runif(K*xcol), K, xcol)
#   q[1,]<-rep(0.01,xcol)
#
#   ## compute loglikelihood
#   for(i in 1:xcol) {
#     f0<-type[[i]][[1]]
#     f0param<-type[[i]][[2]]
#     f1<-type[[i]][[3]]
#     f1param<-type[[i]][[4]]
#     loglike0[[i]]<-f0(x[,i],f0param)
#     loglike1[[i]]<-f1(x[,i],f1param)
#   }
#
#   ## EM algorithm to get MLE of p and q
#   condlike<-list()
#   for(i in 1:xcol) {
#     condlike[[i]]<-matrix(0,xrow,K)
#   }
#
#   loglike.old <- -1e10
#   for(i.iter in 1:max.iter) {
#     if((i.iter%%50) == 0) {
#       print(paste("We have run the first ", i.iter, " iterations for K=", K,sep=""))
#       #print(loglike.old)
#     }
#     err<-tol+1
#
#     ## compute posterior cluster membership
#     clustlike<-matrix(0,xrow,K)
#     templike <- matrix(0,xrow,2)
#     for(j in 1:K) {
#       for(i in 1:xcol) {
#         templike[,1]<-log(q[j,i])+loglike1[[i]]
#         templike[,2]<-log(1-q[j,i])+loglike0[[i]]
#         tempmax<-pmax(templike[,1],templike[,2])
#         for(z in 1:2) {
#           templike[,z]<-exp(templike[,z]-tempmax)
#         }
#         tempsum<-templike[,1]+templike[,2]
#         clustlike[,j]<-clustlike[,j]+tempmax+log(tempsum)
#         condlike[[i]][,j]<-templike[,1]/tempsum
#       }
#       clustlike[,j]<-clustlike[,j]+log(p[j])
#     }
#
#     tempmax<-apply(clustlike,1,max)
#     for(j in 1:K) {
#       clustlike[,j]<-exp(clustlike[,j]-tempmax)
#     }
#     tempsum<-apply(clustlike,1,sum)
#
#
#     ## update motif occurrence rate
#     for(j in 1:K) {
#       clustlike[,j]<-clustlike[,j]/tempsum
#     }
#
#     p.new<-(apply(clustlike,2,sum)+1)/(xrow+K)
#
#     ## update motifs
#     q.new<-matrix(0, K, xcol)
#     for(j in 1:K) {
#       clustpsum<-sum(clustlike[,j])
#       for(i in 1:xcol) {
#         q.new[j,i]<-(sum(clustlike[,j]*condlike[[i]][,j])+1)/(clustpsum+2)
#       }
#     }
#
#     ## evaluate convergence
#     err.p<-max(abs(p.new-p)/p)
#     err.q<-max(abs(q.new-q)/q)
#     err<-max(err.p, err.q)
#
#     ## evaluate whether the log.likelihood increases
#     loglike.new<-(sum(tempmax+log(tempsum))+sum(log(p.new))+sum(log(q.new)+log(1-q.new)))/xrow
#
#
#     p<-p.new
#     q<-q.new
#     loglike.old<-loglike.new
#
#     if(err<tol) {
#       break;
#     }
#   }
#   ## compute posterior p
#   clustlike<-matrix(0,xrow,K)
#   for(j in 1:K) {
#     for(i in 1:xcol) {
#       templike[,1]<-log(q[j,i])+loglike1[[i]]
#       templike[,2]<-log(1-q[j,i])+loglike0[[i]]
#       tempmax<-pmax(templike[,1],templike[,2])
#       for(z in 1:2) {
#         templike[,z]<-exp(templike[,z]-tempmax)
#       }
#       tempsum<-templike[,1]+templike[,2]
#       clustlike[,j]<-clustlike[,j]+tempmax+log(tempsum)
#       condlike[[i]][,j]<-templike[,1]/tempsum
#     }
#     clustlike[,j]<-clustlike[,j]+log(p[j])
#   }
#
#   tempmax<-apply(clustlike,1,max)
#   for(j in 1:K) {
#     clustlike[,j]<-exp(clustlike[,j]-tempmax)
#   }
#   tempsum<-apply(clustlike,1,sum)
#   for(j in 1:K) {
#     clustlike[,j]<-clustlike[,j]/tempsum
#   }
#
#   p.post<-matrix(0,xrow,xcol)
#   for(j in 1:K) {
#     for(i in 1:xcol) {
#       p.post[,i]<-p.post[,i]+clustlike[,j]*condlike[[i]][,j]
#     }
#   }
#
#   ## return
#   #calculate back loglikelihood
#   loglike.old<-loglike.old-(sum(log(p))+sum(log(q)+log(1-q)))/xrow
#   loglike.old<-loglike.old*xrow
#   result<-list(p.post=p.post, motif.prior=p, motif.q=q, loglike=loglike.old,
#                clustlike=clustlike)
# }
#
# generatetype<-function(limfitted)
# {
#   jtype<-list()
#   df<-limfitted$g1num+limfitted$g2num-2+limfitted$df0
#   for(j in 1:limfitted$compnum)
#   {
#     jtype[[j]]<-list(f0=modt.f0.loglike, f0.param=df[j], f1=modt.f1.loglike, f1.param=c(df[j],limfitted$g1num[j],limfitted$g2num[j],limfitted$v0[j]))
#   }
#   jtype
# }
# cormotiffit<-function(exprs = NULL,groupid = NULL,compid = NULL,K=1, tol=1e-3,
#                       max.iter=100,BIC=TRUE, custom_fit = NULL)
# {
#   if (!is.null(custom_fit)) {
#     limfitted <- custom_fit
#     exprs <- matrix(nrow = nrow(limfitted$t))
#   } else if (!any(sapply(list(exprs, groupid, compid), is.null))) {
#     limfitted<-limmafit(exprs,groupid,compid)
#   } else {
#     stop("Improper input")
#   }
#   jtype<-generatetype(limfitted)
#   fitresult<-list()
#   for(i in 1:length(K))
#     fitresult[[i]]<-cmfit(limfitted$t,type=jtype,K=K[i],max.iter=max.iter,tol=tol)
#   bic<-rep(0,length(K))
#   aic<-rep(0,length(K))
#   loglike<-rep(0,length(K))
#   for(i in 1:length(K))
#     loglike[i]<-fitresult[[i]]$loglike
#   for(i in 1:length(K))
#     bic[i]<--2*fitresult[[i]]$loglike+(K[i]-1+K[i]*limfitted$compnum)*log(dim(exprs)[1])
#   for(i in 1:length(K))
#     aic[i]<--2*fitresult[[i]]$loglike+2*(K[i]-1+K[i]*limfitted$compnum)
#   if(BIC==TRUE)
#   {
#     bestflag=which(bic==min(bic))
#   }
#   else
#   {
#     bestflag=which(aic==min(aic))
#   }
#   result<-list(bestmotif=fitresult[[bestflag]],bic=cbind(K,bic),
#                aic=cbind(K,aic),loglike=cbind(K,loglike))
#
# }
#
# plotIC<-function(fitted_cormotif)
# {
#   oldpar<-par(mfrow=c(1,2))
#   plot(fitted_cormotif$bic[,1], fitted_cormotif$bic[,2], type="b",xlab="Motif Number", ylab="BIC", main="BIC")
#   plot(fitted_cormotif$aic[,1], fitted_cormotif$aic[,2], type="b",xlab="Motif Number", ylab="AIC", main="AIC")
#   par(oldpar)
# }
# plotMotif<-function(fitted_cormotif,title="")
# {
#   layout(matrix(1:2,ncol=2))
#   u<-1:dim(fitted_cormotif$bestmotif$motif.q)[2]
#   v<-1:dim(fitted_cormotif$bestmotif$motif.q)[1]
#   image(u,v,t(fitted_cormotif$bestmotif$motif.q),
#         col=gray(seq(from=1,to=0,by=-0.1)),xlab="Study",yaxt = "n",
#         ylab="Corr. Motifs",main=paste(title,"pattern",sep=" "))
#   axis(2,at=1:length(v))
#   for(i in 1:(length(u)+1))
#   {
#     abline(v=(i-0.5))
#   }
#   for(i in 1:(length(v)+1))
#   {
#     abline(h=(i-0.5))
#   }
#   Ng=10000
#   if(is.null(fitted_cormotif$bestmotif$p.post)!=TRUE)
#     Ng=nrow(fitted_cormotif$bestmotif$p.post)
#   genecount=floor(fitted_cormotif$bestmotif$motif.p*Ng)
#   NK=nrow(fitted_cormotif$bestmotif$motif.q)
#   plot(0,0.7,pch=".",xlim=c(0,1.2),ylim=c(0.75,NK+0.25),
#        frame.plot=FALSE,axes=FALSE,xlab="No. of genes",ylab="", main=paste(title,"frequency",sep=" "))
#   segments(0,0.7,fitted_cormotif$bestmotif$motif.p[1],0.7)
#   rect(0,1:NK-0.3,fitted_cormotif$bestmotif$motif.p,1:NK+0.3,
#        col="dark grey")
#   mtext(1:NK,at=1:NK,side=2,cex=0.8)
#   text(fitted_cormotif$bestmotif$motif.p+0.15,1:NK,
#        labels=floor(fitted_cormotif$bestmotif$motif.p*Ng))
# }
