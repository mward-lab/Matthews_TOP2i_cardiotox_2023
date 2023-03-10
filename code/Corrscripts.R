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
cormotif_inital <- cormotiffit(exprs = y_TMM_cpm,
                             groupid = groupid,
                             compid = compid_tran,
                             K=5, max.iter = 500)

cormotif_initial <- cormotif_inital
plotIC(cormotif_initial)

plotMotif(cormotif_initial)

gene_prob_tran <- cormotif_initial$bestmotif$p.post
rownames(gene_prob_tran) <- rownames(y_TMM_cpm)
dim(gene_prob_tran)
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
write_csv(as.data.frame(nonresponse_cluster), "data/cormotif_NRset.txt")
response_cluster24h  <- rownames(gene_prob_tran[(gene_prob_tran[,1]>0.5 &
                                                   gene_prob_tran[,2] >0.5 &
                                                   gene_prob_tran[,3] >0.5 &
                                                   gene_prob_tran[,4] >0.5&
                                                   gene_prob_tran[,5] >0.5 &
                                                   gene_prob_tran[,6] >0.5 &
                                                   gene_prob_tran[,7]>0.5 &
                                                   gene_prob_tran[,8]>0.5 &
                                                   gene_prob_tran[,9]>0.5 &
                                                   gene_prob_tran[,10]>0.5),])
#
 length(response_cluster24h)

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
plotIC(cormotif_24h2mot)
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



