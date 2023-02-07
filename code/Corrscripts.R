# Cormotif
library(edgeR)
library(Cormotif)
library(RColorBrewer)


genas(efit2, coef=c(1,6), subset="Fpval", plot=TRUE, alpha=0.4)
genas(efit2, coef=c(2,7), subset="Fpval", plot=TRUE, alpha=0.4)
genas(efit2, coef=c(3,8), subset="Fpval", plot=TRUE, alpha=0.4)
genas(efit2, coef=c(4,9), subset="Fpval", plot=TRUE, alpha=0.4)
genas(efit2, coef=c(5,10), subset="Fpval", plot=TRUE, alpha=0.4)
genas(efit2, coef=c(1,2), subset="Fpval", plot=TRUE, alpha=0.4)
genas(efit2, coef=c(1,3), subset="Fpval", plot=TRUE, alpha=0.4)
genas(efit2, coef=c(1,4), subset="Fpval", plot=TRUE, alpha=0.4)
genas(efit2, coef=c(1,5), subset="Fpval", plot=TRUE, alpha=0.4)
genas(efit2, coef=c(6,8), subset="Fpval", plot=TRUE, alpha=0.4)
genas(efit2, coef=c(7,9), subset="Fpval", plot=TRUE, alpha=0.4)
genas(efit2, coef=c(8,9), subset="Fpval", plot=TRUE, alpha=0.4)
genas(efit2, coef=c(6,10), subset="Fpval", plot=TRUE, alpha=0.4)




group_fac <- group1
groupid <- as.numeric(group_fac)

compid_tran <- data.frame(c1= c(1,2,3,4,5,7,8,9,10,11), c2 = c( 6,6,6,6,6,12,12,12,12,12))

y_TMM_cpm <- cpm(x, log = TRUE)
colnames(y_TMM_cpm) <- label
y_TMM_cpm
set.seed(12345)
cormotif_tran12 <- cormotiffit(exprs = y_TMM_cpm,
                             groupid = groupid,
                             compid = compid_tran,
                             K=1:12, max.iter = 500)




plotIC(cormotif_tran12)
plotMotif(cormotif_tran12)

head(cormotif_tran$bestmotif$p.post)

gene_prob_tran <- cormotif_tran$bestmotif$p.post
rownames(gene_prob_tran) <- rownames(y_TMM_cpm)
dim(gene_prob_tran)
nonresponse_cluster  <- rownames(gene_prob_tran[(gene_prob_tran[,1] <0.5 & gene_prob_tran[,2] <0.5 & gene_prob_tran[,3] <0.5 & gene_prob_tran[,4] <0.5& gene_prob_tran[,5] <0.5 & gene_prob_tran[,6] <0.5 & gene_prob_tran[,7] <0.5 & gene_prob_tran[,8] <0.5 & gene_prob_tran[,9] <0.5 & gene_prob_tran[,10] <0.5),])

length(nonresponse_cluster)

response_cluster  <- rownames(gene_prob_tran[(gene_prob_tran[,1] >0.5 & gene_prob_tran[,2] >0.5 & gene_prob_tran[,3] >0.5 & gene_prob_tran[,4] >0.5& gene_prob_tran[,5] >0.5 & gene_prob_tran[,6] >0.5 & gene_prob_tran[,7] >0.5 & gene_prob_tran[,8] >0.5 & gene_prob_tran[,9] >0.5 & gene_prob_tran[,10] >0.5),])

length(response_cluster)

