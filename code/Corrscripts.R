# Cormotif
library(edgeR)
library(Cormotif)
library(RColorBrewer)
## read in count file##
design <- read.csv("data/data_outline.txt", row.names = 1)
x_counts <- read.csv("data/norm_counts.csv",row.names = 1)




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

compid_tran <- data.frame(c1= c(1,2,3,4,5,7,8,9,10,11), c2 = c( 6,6,6,6,6,12,12,12,12,12))

y_TMM_cpm <- cpm(x_counts, log = TRUE)
colnames(y_TMM_cpm) <- label
y_TMM_cpm
set.seed(12345)
cormotif_initial <- cormotiffit(exprs = y_TMM_cpm,
                             groupid = groupid,
                             compid = compid_tran,
                             K=5, max.iter = 500)

#saveRDS(cormotif_initial,"data/cormotif_initialK5.RDS")##saved so Ido not have to run everytime
cormotif_initial <- readRDS("data/cormotif_initialK5.RDS")

plotIC(cormotif_initial)

plotMotif(cormotif_initial)

gene_prob_tran <- cormotif_initial$bestmotif$p.post
rownames(gene_prob_tran) <- rownames(y_TMM_cpm)
dim(gene_prob_tran)



#saveRDS(gene_prob_tran,"data/gene_probabilityk5.RDS")

nonresponse_cluster  <- rownames(gene_prob_tran[(gene_prob_tran[,1] <0.5 &
                                                   gene_prob_tran[,2] <0.5 &
                                                   gene_prob_tran[,3] <0.5 &
                                                   gene_prob_tran[,4] <0.5&
                                                   gene_prob_tran[,5] <0.5 &
                                                   gene_prob_tran[,6] <0.5 &
                                                   gene_prob_tran[,7] <0.5 &
                                                   gene_prob_tran[,8] <0.5 &
                                                   gene_prob_tran[,9] <0.5 &
                                                   gene_prob_tran[,10] <0.5),])

length((nonresponse_cluster))
sample(nonresponse_cluster,4)
#"55686"     "100128682" "11190"     "58472"

write_csv(as.data.frame(nonresponse_cluster), "data/cormotif_NRset.txt")

# TI_respset --------------------------------------------------------------
TI_respset  <- rownames(gene_prob_tran[(gene_prob_tran[,1]>0.5 |
                                                   gene_prob_tran[,2] >0.5 |
                                                   gene_prob_tran[,3] >0.5 |
                                                   gene_prob_tran[,4] >0.5|
                                                   # gene_prob_tran[,5] >0.5 |
                                                   gene_prob_tran[,6] >0.5 |
                                                   gene_prob_tran[,7]>0.5 |
                                                   gene_prob_tran[,8]>0.5 |
                                                   gene_prob_tran[,9]>0.5),])
                                                   #gene_prob_tran[,10]>0.5)
#
length(unique(TI_respset))

#5977

write_csv(as.data.frame(TI_respset), "data/cormotif_TI_respset.txt")

TI_respint  <- rownames(gene_prob_tran[(gene_prob_tran[,1]>0.5 &
                                         gene_prob_tran[,2] >0.5 &
                                         gene_prob_tran[,3] >0.5 &
                                         gene_prob_tran[,4] >0.5&
                                         # gene_prob_tran[,5] >0.5 &
                                         gene_prob_tran[,6] >0.5 &
                                         gene_prob_tran[,7]>0.5 &
                                         gene_prob_tran[,8]>0.5 &
                                         gene_prob_tran[,9]>0.5),])

 length(unique(TI_respint))

#95

 write_csv(as.data.frame(TI_respint), "data/cormotif_TI_respint.txt")

 ##eyeball of

# ER_resp set ------------------------------------------------------

 ER_respset  <- rownames(gene_prob_tran[(gene_prob_tran[,1]>0.5 |
                                           gene_prob_tran[,2] >0.5 |
                                           gene_prob_tran[,3] >0.5 |
                                           gene_prob_tran[,4] >0.5),])
                                           # gene_prob_tran[,5] >0.5 |
                                           # gene_prob_tran[,6] >0.5 |
                                           # gene_prob_tran[,7]>0.5 |
                                           # gene_prob_tran[,8]>0.5 |
                                           # gene_prob_tran[,9]>0.5),])
 #gene_prob_tran[,10]>0.5)
 #
 length(unique(ER_respset))

 #916

 write_csv(as.data.frame(ER_respset), "data/cormotif_ER_respset.txt")

 ER_respint  <- rownames(gene_prob_tran[(gene_prob_tran[,1]>0.5 &
                                           gene_prob_tran[,2] >0.5 &
                                           gene_prob_tran[,3] >0.5 &
                                           gene_prob_tran[,4] >0.5),])
                                           # gene_prob_tran[,5] >0.5 &
                                           # gene_prob_tran[,6] >0.5 &
                                           # gene_prob_tran[,7]>0.5 &
                                           # gene_prob_tran[,8]>0.5 &
                                           # gene_prob_tran[,9]>0.5),])
                                           #
 length(unique(ER_respint))

 #898

 write_csv(as.data.frame(ER_respint), "data/cormotif_ER_respint.txt")



# LR_respset --------------------------------------------------------------

 LR_respset  <- rownames(gene_prob_tran[(#gene_prob_tran[,1]>0.5 |
                                           # gene_prob_tran[,2] >0.5 |
                                           # gene_prob_tran[,3] >0.5 |
                                           # gene_prob_tran[,4] >0.5|
                                           # gene_prob_tran[,5] >0.5 |
                                           gene_prob_tran[,6] >0.5 |
                                           gene_prob_tran[,7]>0.5 |
                                           gene_prob_tran[,8]>0.5 |
                                           gene_prob_tran[,9]>0.5),])
 #gene_prob_tran[,10]>0.5)
 #
 length(unique(LR_respset))

 #5515

 write_csv(as.data.frame(LR_respset), "data/cormotif_LR_respset.txt")

 LR_respint  <- rownames(gene_prob_tran[(#gene_prob_tran[,1]>0.5 &
                                           # gene_prob_tran[,2] >0.5 &
                                           # gene_prob_tran[,3] >0.5 &
                                           # gene_prob_tran[,4] >0.5&
                                           # gene_prob_tran[,5] >0.5 &
                                           gene_prob_tran[,6] >0.5 &
                                           gene_prob_tran[,7]>0.5 &
                                           gene_prob_tran[,8]>0.5 &
                                           gene_prob_tran[,9]>0.5),])

 length(unique(LR_respint))

 #1505

 write_csv(as.data.frame(LR_respint), "data/cormotif_LR_respint.txt")



# list summary ------------------------------------------------------------

nonresponse_cluster
 #motif is 8705  set of gene with less than 50% chance of being DE

 TI_respint
 #95
 #intersection of all TOP2Bi genes that have greater than 50% chance of being DE
 TI_respset
 #union of all TOP2Bi genes that have greater than 50% chance of being DE
 #5977

 ER_respint
 #898
 #intersection of all early TOP2Bi genes that have greater than 50% chance of being DE
 ER_respset
 #916
 #union of all early TOP2Bi genes that have greater than 50% chance of being DE



 LR_respint
 #1505
 #intersection of all TOP2Bi genes that have greater than 50% chance of being DE
 LR_respset
 #5515
 #union of all TOP2Bi genes that have greater than 50% chance of being DE


 motif3_TI <- intersect(ER_respset,LR_respset)
 #454

 setdiff(ER_respint,LR_respint)
 #803 genes


 ER_onlyset <- setdiff(ER_respset,LR_respset)
 #462 genes
 motif5_ER <- setdiff(ER_respset,LR_respset)

 motif4_LR <- setdiff(LR_respset,ER_respset)
 #5061
 motif1_NR <- nonresponse_cluster

 DEG_cormotif <- list(motif1_NR,motif3_TI,motif4_LR,motif5_ER)

 names(DEG_cormotif) <- c('motif1_NR','motif3_TI','motif4_LR','motif5_ER')

 saveRDS(DEG_cormotif, "data/DEG_cormotif.RDS")











# old ideas-not in use as of march 20 --------------------------------------




#generesponse tran contains all probabilities

 #use cormotif_initial$bestmotif$motif.q for  probabilities to subset the columns by

q5geneset <- rownames(gene_prob_tran[(gene_prob_tran[,1]> q5[1]|
                            gene_prob_tran[,2] >q5[2]|
                            gene_prob_tran[,3] >q5[3]|
                            gene_prob_tran[,4] >q5[4]),])
                            gene_prob_tran[,5] >q5[5]|
                            gene_prob_tran[,6] >q5[6]|
                            gene_prob_tran[,7]>q5[7]|
                            gene_prob_tran[,8]>q5[8]|
                            gene_prob_tran[,9]>q5[9]|
                            gene_prob_tran[,10]>q5[10]),])

setdiff(q3geneset,q4geneset,)
 list(list_ofrank[which(list_ofrank == 1,arr.ind = TRUE)[,1],])
 row.names(df1[which(df1==17,arr.ind=T)[,1],] )
q1 <- c(0.0003168003, 0.0002657065, 0.0002656879, 0.0002307255, 0.0016080234, 0.0006230686, 0.0002039778,
  0.0002096809, 0.0056795430, 0.0018013869)

 q
 q2 <-  c(0.5310902, 0.4860155, 0.5087463, 0.4850026, 0.4970637, 0.5364329, 0.4777606, 0.4919183, 0.4898462, 0.4964280)
 q3 <- c(0.99633099, 0.98257483, 0.99611512, 0.98937588, 0.08677841, 0.99733673, 0.99080481, 0.99554929, 0.27230843,
 0.04276068)
 q4 <- c(0.0006492446, 0.0004679709, 0.0004171111, 0.0004381761, 0.0028035561, 0.9997297930, 0.9993534905,
0.9993993665, 0.3425304786, 0.0038528269)
 q5 <- c(0.997035453, 0.975360840, 0.996590540, 0.984445187, 0.050353546, 0.010260571, 0.003088745, 0.003893476,
0.017571860, 0.033820708)
 cormotif_initial$bestmotif$motif.q[5,]
 list_ofrank <- generank(cormotif_initial$bestmotif$p.post)

 ### none work so far do eyeballing
 # set1- no response_cluster24h 11190
 # set3- all top2bi 1275442
 # set4-late responsetop2bi only 126820
 # set5-early response top2bi only 27245
 #
 # geneexpressionsets <- cbind(sets=c('set1', 'set3','set4','set5'), ENTREZID = c(11190,1275442,126820,27245))


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
plotMotif(cormotif_24h2mot)
 gene_prob_tran24h <- cormotif_24h2mot$bestmotif$p.post
rownames(gene_prob_tran24h) <- rownames(y_TMM_cpm)
dim(gene_prob_tran24h)
nonresponse_cluster24h <- rownames(gene_prob_tran24h[(gene_prob_tran24h[,1] <0.5 &
                                                        gene_prob_tran24h[,2] <0.5 &
                                                        gene_prob_tran24h[,3] <0.5 &
                                                        gene_prob_tran24h[,4] <0.5&
                                                        gene_prob_tran24h[,5] <0.5),])
length(nonresponse_cluster24h)
response_cluster24h <- rownames(gene_prob_tran24h[(gene_prob_tran24h[,1] > 0.5 &
                                                     gene_prob_tran24h[,2] > 0.5 &
                                                     gene_prob_tran24h[,3] > 0.5 &
                                                     gene_prob_tran24h[,4] > 0.5
                                                    ),])

ACresponse_cluster24h <- rownames(gene_prob_tran24h[(gene_prob_tran24h[,1] >0.5 & gene_prob_tran24h[,2] >0.5 & gene_prob_tran24h[,3] >0.5),])
length(response_cluster24h)
length(ACresponse_cluster24h)
write.csv(response_cluster24h,"data/Top2biresp_cluster24h.csv")
response_cluster24h <- read_csv("data/Top2biresp_cluster24h.csv")
length(intersect(allresponse_cluster,response_cluster24h[[2]]))


write.csv(ACresponse_cluster24h,"data/resp_cluster24h.csv")
write.csv(nonresponse_cluster24h,"data/nonresponse_cluster24h.csv")
# AC3 hour ----------------------------------------------------------------

compid_AC <- data.frame(c1= c(1,2,3), c2 = c( 6,6,6))

y_TMM_cpm <- cpm(x_counts, log = TRUE)
colnames(y_TMM_cpm) <- label
y_TMM_cpm
set.seed(12345)
cormotif_AC3 <- cormotiffit(exprs = y_TMM_cpm,
                               groupid = groupid,
                               compid = compid_AC,
                               K=1:5, max.iter = 500)
plotIC(cormotif_AC3)
plotMotif(cormotif_AC3)


# 24h AC only -------------------------------------------------------------



compid_AC24 <- data.frame(c1= c(7,8,9,10), c2 = c( 12,12,12,12))

y_TMM_cpm <- cpm(x_counts, log = TRUE)
colnames(y_TMM_cpm) <- label
y_TMM_cpm
set.seed(12345)
cormotif_AC24 <- cormotiffit(exprs = y_TMM_cpm,
                               groupid = groupid,
                               compid = compid_AC24,
                               K=1:8, max.iter = 500)
plotIC(cormotif_AC24)
plotMotif(cormotif_AC24)



