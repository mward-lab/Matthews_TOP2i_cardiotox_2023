##BioMart and gene conversion from Dave Tang's Blog'
#https://davetang.org/muse/2013/11/25/thoughts-converting-gene-identifiers/

#install if necessary
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")

library("biomaRt")
library("tidyverse")
library(VennDiagram)
library(kableExtra)
library(scales)
library(ggVennDiagram)
#use the ensembl mart and the human dataset
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
#entrez <- useMart("entrezgene",dataset = "hsapiens_gene_entrez")
#create a filter for all assembled human chromosomes
my_chr <- c(1:22, 'M', 'X', 'Y')


#listAttributes shows all attributes
attributes <- listAttributes(ensembl)
attributes[c(62,63),]
#find entrez attribute name
grep(pattern="HGNC", x=attributes$description, ignore.case=T)
# [1] 58 77 78 79
attributes[c(58, 77:79),]
# name                                 description         page
# 58  entrezgene_trans_name               EntrezGene transcript name ID feature_page
# 77 entrezgene_description NCBI gene (formerly Entrezgene) description feature_page
# 78   entrezgene_accession   NCBI gene (formerly Entrezgene) accession feature_page
# 79          entrezgene_id          NCBI gene (formerly Entrezgene) ID feature_page

#find refseq attribute name
grep(pattern="refseq", x=attributes$description, ignore.case=T)
#[1] 22 23 84 85 86 87 88 89

attributes[c(22,23,84:89),]## what I really want is 84:89
# name                                  description         page
# 22        transcript_mane_select        RefSeq match transcript (MANE Select) feature_page
# 23 transcript_mane_plus_clinical RefSeq match transcript (MANE Plus Clinical) feature_page
# 84                   refseq_mrna                               RefSeq mRNA ID feature_page
# 85         refseq_mrna_predicted                     RefSeq mRNA predicted ID feature_page
# 86                  refseq_ncrna                              RefSeq ncRNA ID feature_page
# 87        refseq_ncrna_predicted                    RefSeq ncRNA predicted ID feature_page
# 88                refseq_peptide                            RefSeq peptide ID feature_page
# 89      refseq_peptide_predicted                  RefSeq peptide predicted ID feature_page
#find ucsc attribute name
#ucsc
grep(pattern="ucsc", x=attributes$description, ignore.case=T)
attributes[94,]
#   name    description         page
# 94 ucsc UCSC Stable ID feature_page

#find Ensembl gene name
attributes[grep(pattern="ensembl",
                     x=attributes,
                     ignore.case=T),]


##has changed to:
# name            description         page
# 1         ensembl_gene_id         Gene stable ID feature_page
# 2 ensembl_gene_id_version Gene stable ID version feature_page

my_refseq_mrna <- getBM(attributes = c('refseq_mrna','refseq_ncrna'),
                        filters = 'chromosome_name',
                        values = my_chr,
                        mart = ensembl)


my_ensembl_gene_id <- getBM(attributes = 'ensembl_gene_id',
                            filters = 'chromosome_name',
                            values = my_chr,
                            mart = ensembl
)


Knowles_2018.elife.33480.supp5.v2 <- read.delim("~/Ward Lab/Cardiotoxicity/Manuscript/Knowles_2018-elife-33480-supp5-v2/Knowles_2018-elife-33480-supp5-v2")
#  read in supp file
knowles5 <- Knowles_2018.elife.33480.supp5.v2[,1]
length(knowles5)  #376 ENSG numbers

response_cluster24h


resp24list <- response_cluster24h$x   #1489  Entrezid
length((resp24list))

resp24list %in% my_full_ensembl_entrez

length(c(supp5list,resp24list) %in% my_full_ensembl_entrez)

# finding names for gene lists ---------------------------------------------
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
my_attributes <- c('entrezgene_id', 'ensembl_gene_id', 'hgnc_symbol')
knowles5 <- Knowles_2018.elife.33480.supp5.v2[,1]
knowles5 <- getBM(attributes=my_attributes,filters ='ensembl_gene_id',
                 values = knowles5, mart = ensembl)
knowles6 <- getBM(attributes=my_attributes,filters ='ensembl_gene_id',
                  values = supp6list, mart = ensembl)
motif24h <- getBM(attributes=my_attributes,filters ='entrezgene_id',
                             values = resp24list, mart = ensembl)
#
# now to make a Venn diagram
length(unique(motif24h$ensembl_gene_id))

# comp knowl v 24h motif --------------------------------------------------



comp1 <- list(motif24h$ensembl_gene_id, knowles5$ensembl_gene_id)
length(unique(c(motif24h$ensembl_gene_id, knowles5$ensembl_gene_id)))
ggVennDiagram(comp1,
              category.names = c("Top2Bi-24hours","Knowles-supp5"),
              show_intersect = FALSE,
              set_color = "black",
              label = "both",
              label_percent_digit = 1,
              label_size = 3,
              label_alpha = 0,
              label_color = "black",
              edge_lty = "solid", set_size = )+
  scale_x_continuous(expand = expansion(mult = .21))+
  scale_fill_gradient(low = "red", high = "yellow")+
  labs(title = "Comparision Knowles v 24h", caption = "n = 1889")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))
comp2 <- list(motif24h$ensembl_gene_id, knowles6$ensembl_gene_id)
length(unique(c(motif24h$ensembl_gene_id, knowles6$ensembl_gene_id)))

ggVennDiagram(comp2,
              category.names = c("Top2Bi-24hours","Knowles-supp6"),
              show_intersect = FALSE,
              set_color = "black",
              label = "both",
              label_percent_digit = 1,
              label_size = 3,
              label_alpha = 0,
              label_color = "black",
              edge_lty = "solid", set_size = )+
  scale_x_continuous(expand = expansion(mult = .21))+
  scale_fill_gradient(low = "red", high = "yellow")+
  labs(title = "Comparision Knowles v 24h", caption = "n = 1960")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))
# comp knowles v sig diff 24AC --------------------------------------------
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
my_attributes <- c('entrezgene_id', 'ensembl_gene_id', 'hgnc_symbol')
knowles5 <- Knowles_2018.elife.33480.supp5.v2[,1]
knowles5 <- getBM(attributes=my_attributes,filters ='ensembl_gene_id',
                 values = knowles5, mart = ensembl)
in_common24AC <- intersect(sigVDA24$ENTREZID,sigVDX24$ENTREZID)
in_common24AC <- intersect(in_common24AC,sigVEP24$ENTREZID)


AC24hoursig <- (unique(in_common24AC))  ### 5400 total genes from all 24 AC
AConly24hsig <- getBM(attributes=my_attributes,filters ='entrezgene_id',
                  values = AC24hoursig, mart = ensembl)
comp3 <- list(AConly24hsig$ensembl_gene_id, knowles5$ensembl_gene_id)
length(unique(c(AConly24hsig$ensembl_gene_id, knowles5$ensembl_gene_id)))
ggVennDiagram(comp3,
              category.names = c("AC 24 hour only sig","Knowles-supp5"),
              show_intersect = FALSE,
              set_color = "black",
              label = "both",
              label_percent_digit = 1,
              label_size = 3,
              label_alpha = 0,
              label_color = "black",
              edge_lty = "solid", set_size = )+
  scale_x_continuous(expand = expansion(mult = .21))+
  scale_fill_gradient(low = "red", high = "yellow")+
  labs(title = "Comparision Knowles-5 v 24h sig AC", caption = "n = 5958")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

intersect(AConly24hsig$ensembl_gene_id,knowles5$ensembl_gene_id) #154 genes to GO on
intersect(motif24h$ensembl_gene_id,knowles5$ensembl_gene_id) ##47
intersect(AConly24hsig$ensembl_gene_id,knowles6$ensembl_gene_id) #183genes to GO on
intersect(motif24h$ensembl_gene_id,knowles6$ensembl_gene_id) ##57

intersect(AConly24hsig$ensembl_gene_id,motif24h$ensembl_gene_id)  ###1450


NR_cluster24h <-nonresponse_cluster24h$x ### 5400 total genes from all 24 AC
NR_cluster24h <- getBM(attributes=my_attributes,filters ='entrezgene_id',
                      values = NR_cluster24h, mart = ensembl)
intersect(AConly24hsig$ensembl_gene_id,NR_cluster24h$ensembl_gene_id)





comp4 <- list(NR_cluster24h$ensembl_gene_id, knowles5$ensembl_gene_id)
length(unique(c(NR_cluster24h$ensembl_gene_id, knowles5$ensembl_gene_id)))
ggVennDiagram(comp4,
              category.names = c("Motif 24h No-response","Knowles-supp5"),
              show_intersect = FALSE,
              set_color = "black",
              label = "both",
              label_percent_digit = 1,
              label_size = 3,
              label_alpha = 0,
              label_color = "black",
              edge_lty = "solid", set_size = )+
  scale_x_continuous(expand = expansion(mult = .21))+
  scale_fill_gradient(low = "red", high = "yellow")+
  labs(title = "Comparision Knowles-5 v No-response genes", caption = "n = 9934")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

comp5 <-
length(unique(c(AConly24hsig$ensembl_gene_id, motif24h$ensembl_gene_id)))

ggVennDiagram(comp5,
              category.names = c("Motif 24h response","sig 24h DEG AC"),
              show_intersect = FALSE,
              set_color = "black",
              label = "both",
              label_percent_digit = 1,
              label_size = 3,
              label_alpha = 0,
              label_color = "black",
              edge_lty = "solid", set_size = )+
  scale_x_continuous(expand = expansion(mult = .21))+
  scale_fill_gradient(low = "blue", high = "yellow")+
  labs(title = "Comparision Motif 24h response to sigDEG AC drugs", caption = "n = 5855")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

specialset <- intersect(AConly24hsig$ensembl_gene_id, motif24h$ensembl_gene_id)
comp6 <- list(specialset,knowles5$ensembl_gene_id )
length(unique(c(specialset,knowles5$ensembl_gene_id)))

ggVennDiagram(comp6,
              category.names = c("significant and motif overlap","Knowles-supp 5"),
              show_intersect = FALSE,
              set_color = "black",
              label = "both",
              label_percent_digit = 1,
              label_size = 3,
              label_alpha = 0,
              label_color = "black",
              edge_lty = "solid", set_size = )+
  scale_x_continuous(expand = expansion(mult = .21))+
  scale_fill_gradient(low = "blue", high = "yellow")+
  labs(title = "Overlaps of my data sets and knowles supp5", caption = "n = 1784")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))
length(unique(NR_cluster24h$ensembl_gene_id))


sigDox24 <- getBM(attributes=my_attributes,filters ='entrezgene_id',
                                  values = sigVDX24$ENTREZID, mart = ensembl)
comp7 <- list(sigDox24$ensembl_gene_id,knowles5$ensembl_gene_id)
length(unique(c(sigDox24$ensembl_gene_id,knowles5$ensembl_gene_id)))




ggVennDiagram(comp7,
              category.names = c("significant 24hour Doxorubicin","Knowles-supp 5"),
              show_intersect = FALSE,
              set_color = "black",
              label = "both",
              label_percent_digit = 1,
              label_size = 3,
              label_alpha = 0,
              label_color = "black",
              edge_lty = "solid", set_size = )+
  scale_x_continuous(expand = expansion(mult = .21))+
  scale_fill_gradient(low = "lightblue", high = "darkblue")+
  labs(title = "Doxorubicin 24 hour sig genes and knowles supp5", caption = "n = 7373")+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))

# ACresponse cluster ------------------------------------------------------

AC_cluster24h <-ACresponse_cluster24h$x ### 5400 total genes from all 24 AC
AC_cluster24h <- getBM(attributes=my_attributes,filters ='entrezgene_id',
                       values = AC_cluster24h, mart = ensembl)
length(unique(AC_cluster24h$ensembl_gene_id))#5799
intersect(AC_cluster24h$ensembl_gene_id, NR_cluster24h$ensembl_gene_id)#[7]
intersect(AC_cluster24h$ensembl_gene_id, knowles5$ensembl_gene_id)#154
length(intersect(AC_cluster24h$ensembl_gene_id, motif24h$ensembl_gene_id)) #1569
length(intersect(sigDox24$ensembl_gene_id, AC_cluster24h$ensembl_gene_id))#5472
length(intersect(sigDox24$ensembl_gene_id, knowles5$ensembl_gene_id))#185
length(intersect(sigDox24$ensembl_gene_id, motif24h$ensembl_gene_id))#1513
length(intersect(AConly24hsig$ensembl_gene_id, AC_cluster24h$ensembl_gene_id))#5106

length(intersect(sigDox24$ensembl_gene_id, NR_cluster24h$ensembl_gene_id))#1710
length(intersect(motif24h$ensembl_gene_id, NR_cluster24h$ensembl_gene_id))# 2
length(intersect(AConly24hsig$ensembl_gene_id, NR_cluster24h$ensembl_gene_id))#634
length(intersect(AC_cluster24h$ensembl_gene_id, knowles5$ensembl_gene_id))# 154
length(intersect(AC_cluster24h$ensembl_gene_id, motif24h$ensembl_gene_id))#1569

length(intersect(AConly24hsig$ensembl_gene_id, knowles5$ensembl_gene_id))#154

# supp6 -------------------------------------------------------------------

length(intersect(knowles6$ensembl_gene_id, NR_cluster24h$ensembl_gene_id))#260
length(intersect(knowles6$ensembl_gene_id, sigDox24$ensembl_gene_id))#216
length(intersect(knowles6$ensembl_gene_id, AC_cluster24h$ensembl_gene_id))#181
length(intersect(knowles6$ensembl_gene_id, motif24h$ensembl_gene_id))
length(intersect(knowles6$ensembl_gene_id, AConly24hsig$ensembl_gene_id))#183
length(intersect(knowles6$entrezgene_id, NR_cluster24h$entrezgene_id))



# pairwise ----------------------------------------------------------------
# First make a .GMT of the Knowles supp 5 and 6 file.   This would be to pivot wider the 377 and 447 genes reported in file:
knowleslist <- list(knowles5$ensembl_gene_id,knowles6$ensembl_gene_id)

writeGMT(knowleslist, "data/knowlesGMT.GMT")##get it in the write form, but missing two columns added by hand
#response genes:  sigvDX24h  unique ensemble:  7182
#total response genes:

#read in file:

# knowles56 <- read_delim("data/knowles56.GMT",
#                         delim = "\t", escape_double = FALSE,
#                         col_names = FALSE, trim_ws = TRUE)
# library(gprofiler2)
# sDXgostres <- gost(query = sigDox24$ensembl_gene_id,  organism = "hsapiens",
#                                                      ordered_query = FALSE,
#                                                      domain_scope = "custom",
#                                                      measure_underrepresentation = FALSE,
#                                                      evcodes = FALSE,
#                                                      user_threshold = 0.05,
#                                                      correction_method = "fdr",
#                                                      custom_bg = totalgenes$ensembl_gene_id,
#                                                      sources="data/knowles56.GMT")

length(intersect(totalgenes$ensembl_gene_id, knowles5$ensembl_gene_id))#375
length(unique(c(totalgenes$ensembl_gene_id, knowles5$ensembl_gene_id)))

length(sigDox24$ensembl_gene_id %in% !totalgenes$ensembl_gene_id )
write_delim("data/totalgenes.txt",delim = "\t", escape_double = FALSE,
                                  col_names = FALSE, trim_ws = TRUE)
library(readr)

write_tsv(as.data.frame(totalgenes$ensembl_gene_id), file ="data/ensgtotal.txt")



# antijoin: ---------------------------------------------------------------

NRsigDX <- totalgenes %>% anti_join(sigDox24, by = "ensembl_gene_id")
length(intersect(NRsigDX$ensembl_gene_id, knowles5$ensembl_gene_id))
length(unique(NRsigDX$ensembl_gene_id))
NRAConly24hsig <- totalgenes %>% anti_join(AConly24hsig, by = "ensembl_gene_id")
length(unique(NRAConly24hsig$ensembl_gene_id))
length(intersect(NRAConly24hsig$ensembl_gene_id, knowles5$ensembl_gene_id))
length(intersect(AConly24hsig$ensembl_gene_id, knowles5$ensembl_gene_id))




# Keep counting: ----------------------------------------------------------

knowles5
NoResp <- getBM(attributes=my_attributes,filters ='entrezgene_id',
                       values = NoResp$ENTREZID, mart = ensembl)
top2bi <- getBM(attributes=my_attributes,filters ='entrezgene_id',
                values = top2bi$x, mart = ensembl)
Dxresp <- getBM(attributes=my_attributes,filters ='entrezgene_id',
      values = Dxresp, mart = ensembl)
DDEMresp <- getBM(attributes=my_attributes,filters ='entrezgene_id',
      values = DDEMresp, mart = ensembl)
DDEresp <- getBM(attributes=my_attributes,filters ='entrezgene_id',
                  values = DDEresp, mart = ensembl)

length(intersect(knowles5$entrezgene_id,Dxresp$entrezgene_id))#203****NEW IS 185
  reQTLDX <- intersect(knowles5$entrezgene_id, Dxresp$entrezgene_id)
length(intersect(knowles5$entrezgene_id,DDEresp$entrezgene_id))#123
  reQTLDDEint <- intersect(knowles5$entrezgene_id,DDEresp$entrezgene_id)
length(intersect(knowles5$entrezgene_id,DDEMresp$entrezgene_id))# 31
  reQTLDDEMint <- intersect(knowles5$entrezgene_id,DDEMresp$entrezgene_id)
length(intersect(knowles5$entrezgene_id,NoResp$entrezgene_id))#145
  reQTLNRresp <- intersect(knowles5$entrezgene_id,NoResp$entrezgene_id)


daepint <- intersect(sigVDA24$ENTREZID, sigVEP24$ENTREZID)
daepuni <- union(sigVDA24$ENTREZID, sigVEP24$ENTREZID)
DaEprespint <- getBM(attributes=my_attributes,filters ='entrezgene_id',
                values = daepint, mart = ensembl)
DaEprespuni <- getBM(attributes=my_attributes,filters ='entrezgene_id',
                     values = daepuni, mart = ensembl)

AConlyresp <- union(Dxresp$entrezgene_id, DaEprespuni$entrezgene_id)

reQTLDDEunion <- intersect(AConlyresp, knowles5$entrezgene_id)
write_tsv(as.data.frame(reQTLDDEint), file ="data/DDE_reQTL.txt")
reQTLDX <- dx24intknowles

setdiff(reQTLDX, reQTLDDEintersect)
intersect(setdiff(reQTLDX, reQTLDDEintersect),DDEMresp$entrezgene_id)
reQTL_DaEpunion <- intersect(knowles5$entrezgene_id, DaEprespuni$entrezgene_id)
reQTL_DaEpint <-  intersect(knowles5$entrezgene_id, DaEprespint$entrezgene_id)


#reQTLDDEunion
length(intersect(reQTLDDEunion,reQTLDDEunion))
length(intersect(reQTLDDEunion,reQTLDDEint))
length(intersect(reQTLDDEunion,reQTLDX))
length(intersect(reQTLDDEunion,reQTL_DaEpunion))
length(intersect(reQTLDDEunion,reQTL_DaEpint))


#reQTLDDEint
length(intersect(reQTLDDEint,reQTLDDEunion))
length(intersect(reQTLDDEint,reQTLDDEint))
length(intersect(reQTLDDEint,reQTLDX))
length(intersect(reQTLDDEint,reQTL_DaEpunion))
length(intersect(reQTLDDEint,reQTL_DaEpint))

#reQTLDX
length(intersect(reQTLDX,reQTLDDEunion))
length(intersect(reQTLDX,reQTLDDEint))
length(intersect(reQTLDX,reQTLDX))
length(intersect(reQTLDX,reQTL_DaEpunion))
length(intersect(reQTLDX,reQTL_DaEpint))

#reQTL_DaEpunion
length(intersect(reQTL_DaEpunion,reQTLDDEunion))
length(intersect(reQTL_DaEpunion,reQTLDDEint))
length(intersect(reQTL_DaEpunion,reQTLDX))
length(intersect(reQTL_DaEpunion,reQTL_DaEpunion))
length(intersect(reQTL_DaEpunion,reQTL_DaEpint))

#reQTL_DaEpint
length(intersect(reQTL_DaEpint,reQTLDDEunion))
length(intersect(reQTL_DaEpint,reQTLDDEint))
length(intersect(reQTL_DaEpint,reQTLDX))
length(intersect(reQTL_DaEpint,reQTL_DaEpunion))
length(intersect(reQTL_DaEpint,reQTL_DaEpint))


# setdiff -----------------------------------------------------------------
#reQTLDDEunion
length(setdiff(reQTLDDEunion,reQTLDDEunion))
length(setdiff(reQTLDDEunion,reQTLDDEint))
length(setdiff(reQTLDDEunion,reQTLDX))
length(setdiff(reQTLDDEunion,reQTL_DaEpunion))
length(setdiff(reQTLDDEunion,reQTL_DaEpint))


#reQTLDDEint
length(setdiff(reQTLDDEint,reQTLDDEunion))
length(setdiff(reQTLDDEint,reQTLDDEint))
length(setdiff(reQTLDDEint,reQTLDX))
length(setdiff(reQTLDDEint,reQTL_DaEpunion))
length(setdiff(reQTLDDEint,reQTL_DaEpint))

#reQTLDX
length(setdiff(reQTLDX,reQTLDDEunion))
length(setdiff(reQTLDX,reQTLDDEint))
length(setdiff(reQTLDX,reQTLDX))
length(setdiff(reQTLDX,reQTL_DaEpunion))
DXspec_reQTLgenes<- setdiff(reQTLDX,reQTL_DaEpunion)
#write.delim(DXspec_reQTLgenes,"data/Dx_reQTL_specific.txt")
write_tsv(as.data.frame(DXspec_reQTLgenes), file ="data/Dx_reQTL_specific.txt")


length(setdiff(reQTLDX,reQTL_DaEpint))

#reQTL_DaEpunion
length(setdiff(reQTL_DaEpunion,reQTLDDEunion))
length(setdiff(reQTL_DaEpunion,reQTLDDEint))
length(setdiff(reQTL_DaEpunion,reQTLDX))
length(setdiff(reQTL_DaEpunion,reQTL_DaEpunion))
length(setdiff(reQTL_DaEpunion,reQTL_DaEpint))

#reQTL_DaEpint
length(setdiff(reQTL_DaEpint,reQTLDDEunion))
length(setdiff(reQTL_DaEpint,reQTLDDEint))
length(setdiff(reQTL_DaEpint,reQTLDX))
length(setdiff(reQTL_DaEpint,reQTL_DaEpunion))
length(setdiff(reQTL_DaEpint,reQTL_DaEpint))





# Only dox gene set table 2-16-2023 -------------------------------------------------------

DXsprespon <- list24totvenn$..values..[[14]]
ONLYDXgenels <- getBM(attributes=my_attributes,filters ='entrezgene_id',
                values = DXsprespon, mart = ensembl)
DXspec_reQTL <- unique(ONLYDXgenels$entrezgene_id)

length(intersect(DXspec_reQTL,reQTLDDEunion))


length(intersect(DXspec_reQTL,reQTLDDEint))
length(intersect(DXspec_reQTL,reQTLDX))
length(intersect(DXspec_reQTL,reQTL_DaEpunion))
length(intersect(Dxrespon,reQTL_DaEpint))
sum(duplicated(ONLYDXgenels$entrezgene_id))
#[1] 18

length(intersect(DXsprespon,knowles5$entrezgene_id))
length(setdiff(knowles5$entrezgene_id,DXsprespon))
DXsprespon
length(unique(ONLYDXgenels))

thing <- setdiff(knowles5$entrezgene_id,DXsprespon)
thing2 <- setdiff(knowles5$entrezgene_id, reQTLDDEMint)
length(intersect (thing, thing2))
setdiff(thing, thing2)
