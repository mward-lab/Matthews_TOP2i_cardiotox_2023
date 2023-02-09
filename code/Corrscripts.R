# Cormotif
library(edgeR)
library(Cormotif)
library(RColorBrewer)
## read in count file##
design <- read.csv("data/data_outline.txt", row.names = 1)
x_counts <- read.csv("data/norm_counts.csv",row.names = 1)

# genas(efit2, coef=c(1,6), subset="Fpval", plot=TRUE, alpha=0.4)
#
# genas(efit2, coef=c(2,7), subset="Fpval", plot=TRUE, alpha=0.4)
# genas(efit2, coef=c(3,8), subset="Fpval", plot=TRUE, alpha=0.4)
# genas(efit2, coef=c(4,9), subset="Fpval", plot=TRUE, alpha=0.4)
# genas(efit2, coef=c(5,10), subset="Fpval", plot=TRUE, alpha=0.4)
# genas(efit2, coef=c(1,2), subset="Fpval", plot=TRUE, alpha=0.4)
# genas(efit2, coef=c(1,3), subset="Fpval", plot=TRUE, alpha=0.4)
# genas(efit2, coef=c(1,4), subset="Fpval", plot=TRUE, alpha=0.4)
# genas(efit2, coef=c(1,5), subset="Fpval", plot=TRUE, alpha=0.4)
# genas(efit2, coef=c(6,8), subset="Fpval", plot=TRUE, alpha=0.4)
# genas(efit2, coef=c(7,9), subset="Fpval", plot=TRUE, alpha=0.4)
# genas(efit2, coef=c(8,9), subset="Fpval", plot=TRUE, alpha=0.4)
# genas(efit2, coef=c(6,10), subset="Fpval", plot=TRUE, alpha=0.4)

colnames(shortcompmat) <- c("Test","Control")
groupmat <- matrix(rep(c("Daunorubicin.24h","Vehicle.24h"),12),nrow =1, ncol=24 )



# All data ----------------------------------------------------------------


group_fac <- group1
groupid <- as.numeric(group_fac)

compid_tran <- data.frame(c1= c(1,2,3,4,5,7,8,9,10,11), c2 = c( 6,6,6,6,6,12,12,12,12,12))

y_TMM_cpm <- cpm(x_counts, log = TRUE)
colnames(y_TMM_cpm) <- label
y_TMM_cpm
set.seed(12345)
cormotif_tran12 <- cormotiffit(exprs = y_TMM_cpm,
                             groupid = groupid,
                             compid = compid_tran,
                             K=5, max.iter = 500)




plotIC(cormotif_tran12)
colnames(cormotif_tran12$bestmotif$motif.q) <- c("3_Daun","3_Dox","3_Epi","3_Mito","3_Tras","24_Daun","24_Dox","24_Epi","24_Mito","24_Tras")
plotMotif(cormotif_tran12)

head(cormotif_tran12$bestmotif$p.post)

gene_prob_tran <- cormotif_tran$bestmotif$p.post
rownames(gene_prob_tran) <- rownames(y_TMM_cpm)
# dim(gene_prob_tran)
# nonresponse_cluster  <- rownames(gene_prob_tran[(gene_prob_tran[,1] <0.5 & gene_prob_tran[,2] <0.5 & gene_prob_tran[,3] <0.5 & gene_prob_tran[,4] <0.5& gene_prob_tran[,5] <0.5 & gene_prob_tran[,6] <0.5 & gene_prob_tran[,7] <0.5 & gene_prob_tran[,8] <0.5 & gene_prob_tran[,9] <0.5 & gene_prob_tran[,10] <0.5),])
#
# length(nonresponse_cluster)
#
# response_cluster  <- rownames(gene_prob_tran[(gene_prob_tran[,1] >0.5 & gene_prob_tran[,2] >0.5 & gene_prob_tran[,3] >0.5 & gene_prob_tran[,4] >0.5& gene_prob_tran[,5] >0.5 & gene_prob_tran[,6] >0.5 & gene_prob_tran[,7] >0.5 & gene_prob_tran[,8] >0.5 & gene_prob_tran[,9] >0.5 & gene_prob_tran[,10] >0.5),])
#
# length(response_cluster)

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


# 24 hour patterns -------------------------------------h3groupid <- factor(rep(1:6,6),levels = c(1,2,3,4,5,6))
compid_24h <- data.frame(c1= c(7,8,9,10,11), c2 = c( 12,12,12,12,12))

y_TMM_cpm <- cpm(x_counts, log = TRUE)


colnames(y_TMM_cpm) <- label
#threehour <- y_TMM_cpm[,c(1:6,13:18,25:30,37:42,49:54,61:66)]
set.seed(12345)
cormotif_24h <- cormotiffit(exprs = y_TMM_cpm,
                           groupid = groupid,
                           compid = compid_24h,
                           K=1:8, max.iter = 500)

plotIC(cormotif_24h)
#x_axis_labels(labels = c("3_Daun","3_Dox","3_Epi","3_Mito","3_Tras","24_Daun","24_Dox","24_Epi","24_Mito","24_Tras"), every_nth = 1, adj=1, srt =90, cex =0.4)
plotMotif(cormotif_24h)


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
                               K=2:8, max.iter = 500)
plotIC(cormotif_AC24)
plotMotif(cormotif_AC24)
